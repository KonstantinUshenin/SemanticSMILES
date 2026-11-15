from __future__ import annotations

from pathlib import Path
from typing import Any, TypedDict, List, Optional, Dict

from langgraph.graph import StateGraph, START, END
from langfuse.langchain import CallbackHandler
from langchain_openai import ChatOpenAI

from pydantic import BaseModel, Field
import os
import sys
import yaml


from dotenv import load_dotenv
load_dotenv()
OPENAI_API_KEY = os.getenv("OPENAI_API_KEY")
OPENAI_BASE_URL = os.getenv("OPENAI_API_BASE")
OPENAI_MODEL = os.getenv("OPENAI_MODEL")

_PATHS_FILE = Path(__file__).resolve().parent.parent / "paths.yaml"
with _PATHS_FILE.open(encoding="utf-8") as _f:
    _PATHS = yaml.safe_load(_f)

WORKSPACE = _PATHS["workspace"]
WORKSPACE_DATA = _PATHS["workspace_data"]
WORKSPACE_SUMMARY = _PATHS["workspace_summary"]

from langfuse import Langfuse

langfuse = Langfuse()
print(langfuse.auth_check())  # should return True

from langfuse.langchain import CallbackHandler
langfuse_handler = CallbackHandler()  # reads env vars

LF_CONFIG: Dict[str, Any] = {
    "callbacks": [langfuse_handler],
}

# Langfuse CallbackHandler uses OpenTelemetry on the graph root; LangGraph parallel ``Send``
# branches end spans in different asyncio tasks, so context detach can raise ValueError and
# print "Failed to detach context". Omit callbacks for fan-out graphs, or use max_concurrency=1
# if you need Langfuse traces for that run.
LF_CONFIG_LANGGRAPH_FANOUT: Dict[str, Any] = {}

LLM_BASIC = ChatOpenAI(
    model=OPENAI_MODEL,
    base_url=OPENAI_BASE_URL,
    api_key=OPENAI_API_KEY,
    temperature=0.1,
)

RETRY = 3






















import asyncio
import threading
import uuid
from dataclasses import asdict, is_dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping

from pydantic import BaseModel


class VersionedYaml:
    """
    Immutable YAML storage by table:
    - path format: db_root/table/YYYYMMDD_NNNN.yaml
    - write-only records (no update/delete)
    """

    def __init__(self, db_root: str | Path):
        self.db_root = Path(db_root)
        self.db_root.mkdir(parents=True, exist_ok=True)
        self._meta_lock = threading.RLock()
        self._table_locks: dict[str, threading.RLock] = {}

    def create(
        self,
        table: str,
        record: Any,
        _meta_surplus: Mapping[str, Any] | None = None,
    ) -> Path:
        table_name = self._validate_table(table)
        lock = self._table_lock(table_name)

        with lock:
            table_dir = self.db_root / table_name
            table_dir.mkdir(parents=True, exist_ok=True)

            now = datetime.now(timezone.utc)
            day = now.strftime("%Y%m%d")
            version = self._next_version_locked(table_dir, day)
            key = f"{day}_{version:04d}"
            path = table_dir / f"{key}.yaml"

            meta = {
                "table": table_name,
                "key": key,
                "day": day,
                "version": version,
                "created_at_utc": now.isoformat(),
            }
            if _meta_surplus:
                for k, v in dict(_meta_surplus).items():
                    if k not in meta:
                        meta[k] = v

            payload = {
                "_meta": meta,
                "_state": self._normalize_record(record),
            }
            self._atomic_dump_yaml(path, payload)
            return path

    async def acreate(
        self,
        table: str,
        record: Any,
        _meta_surplus: Mapping[str, Any] | None = None,
    ) -> Path:
        return await asyncio.to_thread(self.create, table, record, _meta_surplus)

    def read(self, table: str) -> tuple[Any, dict[str, Any]] | None:
        table_name = self._validate_table(table)
        lock = self._table_lock(table_name)

        with lock:
            table_dir = self.db_root / table_name
            latest_path = self._latest_file_locked(table_dir)
            if latest_path is None:
                return None
            payload = self._load_yaml(latest_path)

        record = payload.get("_state", payload.get("state"))
        meta = payload.get("_meta", {})
        if not isinstance(meta, dict):
            meta = {}
        return record, meta

    async def aread(self, table: str) -> tuple[Any, dict[str, Any]] | None:
        return await asyncio.to_thread(self.read, table)

    def _table_lock(self, table: str) -> threading.RLock:
        with self._meta_lock:
            lock = self._table_locks.get(table)
            if lock is None:
                lock = threading.RLock()
                self._table_locks[table] = lock
            return lock

    @staticmethod
    def _validate_table(table: str) -> str:
        table_name = str(table).strip()
        if not table_name:
            raise ValueError("table must be non-empty.")
        return table_name

    def _next_version_locked(self, table_dir: Path, day: str) -> int:
        max_version = -1
        for file_path in table_dir.glob(f"{day}_*.yaml"):
            name = file_path.stem
            parts = name.split("_", 1)
            if len(parts) != 2:
                continue
            p_day, p_ver = parts
            if p_day != day or len(p_ver) != 4 or not p_ver.isdigit():
                continue
            max_version = max(max_version, int(p_ver))
        next_version = max_version + 1
        if next_version > 9999:
            raise ValueError("Too many versions for one day (max 9999).")
        return next_version

    def _latest_file_locked(self, table_dir: Path) -> Path | None:
        if not table_dir.exists():
            return None

        latest_key: tuple[str, int] | None = None
        latest_path: Path | None = None
        for file_path in table_dir.glob("*.yaml"):
            name = file_path.stem
            parts = name.split("_", 1)
            if len(parts) != 2:
                continue
            day, ver = parts
            if len(day) != 8 or not day.isdigit() or len(ver) != 4 or not ver.isdigit():
                continue
            key = (day, int(ver))
            if latest_key is None or key > latest_key:
                latest_key = key
                latest_path = file_path
        return latest_path

    def _normalize_record(self, record: Any) -> Any:
        if isinstance(record, BaseModel):
            return record.model_dump(mode="python")
        if hasattr(record, "model_dump") and callable(record.model_dump):
            return record.model_dump(mode="python")
        if is_dataclass(record):
            return asdict(record)
        if isinstance(record, Mapping):
            return dict(record)
        return record

    @staticmethod
    def _atomic_dump_yaml(path: Path, payload: dict[str, Any]) -> None:
        tmp_path = path.with_suffix(path.suffix + f".tmp-{uuid.uuid4().hex}")
        with tmp_path.open("w", encoding="utf-8") as handle:
            yaml.safe_dump(payload, handle, sort_keys=False, allow_unicode=False)
        tmp_path.replace(path)

    @staticmethod
    def _load_yaml(path: Path) -> dict[str, Any]:
        with path.open("r", encoding="utf-8") as handle:
            payload = yaml.safe_load(handle)
        if payload is None:
            return {}
        if isinstance(payload, dict):
            return payload
        return {"_state": payload}



class _PromptRegistry:
    """Internal registry singleton (do not instantiate manually)."""
    def __init__(self):
        self._registry: Dict[str, PromptEntry] = {}

    def register(
        self,
        name: str,
        template: str,
        output_model: Optional[Type[BaseModel]] = None,
        description: str = "",
        version: str = "v1",
    ):
        if name in self._registry:
            raise ValueError(f"Prompt '{name}' already registered.")
        entry = PromptEntry(name, template, output_model, description, version)
        self._registry[name] = entry
        return entry

    def get(self, name: str) -> PromptEntry:
        try:
            return self._registry[name]
        except KeyError:
            raise KeyError(f"Prompt '{name}' not found in registry.")

    def all(self) -> Dict[str, PromptEntry]:
        return dict(self._registry)


# ✅ Guarantee a single global instance
if "chemdataset.llm.PROMPT_REGISTRY" not in sys.modules:
    sys.modules["chemdataset.llm.PROMPT_REGISTRY"] = _PromptRegistry()

# Expose the singleton directly
PromptRegistry: _PromptRegistry = sys.modules["chemdataset.llm.PROMPT_REGISTRY"]
