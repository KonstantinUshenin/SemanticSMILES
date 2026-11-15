# Create Dataset

Jupyter notebooks that turn a folder of source chemistry datasets into curated SMILES tables, semantic tags, and summary descriptions (the SemanticSMILES construction pipeline).

Each source dataset lives under `workspace/<DATASET_KEY>/` (README, CSV files, etc.). Intermediate and final products are written to `workspace_DATA/` and `workspace_SUMMARY/`.

## Layout

```
code/
├── paths.yaml          # workspace path config
├── .env.example        # LLM / Langfuse credentials template
├── README.md
└── notebooks/
    ├── core.py         # LLM client, Langfuse, VersionedYaml, path loading
    ├── processing.py   # SMILES parsing / BasicAnalysis
    ├── analysis.py     # RDKit-based molecular feature helpers
    ├── structure_utils.py
    ├── utils.py
    ├── row.py
    ├── llm_check.py
    ├── step01_detect_files.ipynb
    ├── step02_detect_tags.ipynb
    ├── step03_detect_SMILES.ipynb
    ├── step04_process_SMILES.ipynb
    └── step05_process_SMILES.ipynb
```

## Prerequisites

- Python 3.11+ (recommended)
- Jupyter / VS Code / Cursor notebook runner
- OpenAI-compatible LLM endpoint (API key + base URL)
- Optional: Langfuse for tracing (`LANGFUSE_*` env vars)
- RDKit (steps 04–05), pandas, pydantic, PyYAML, python-dotenv
- LangChain / LangGraph / `langchain-openai` / langfuse (steps 01–03)

Install packages into your environment as needed, for example:

```bash
pip install jupyter pandas pydantic pyyaml python-dotenv rdkit \
  langchain langgraph langchain-openai langfuse
```

## Configuration

### 1. Paths — [`paths.yaml`](paths.yaml)

Edit this file before running notebooks. Values are loaded by [`notebooks/core.py`](notebooks/core.py) and exported as `WORKSPACE`, `WORKSPACE_DATA`, and `WORKSPACE_SUMMARY`.

```yaml
workspace: /path/to/workspace
workspace_data: /path/to/workspace_DATA
workspace_summary: /path/to/workspace_SUMMARY
```

| Key | Role |
|-----|------|
| `workspace` | Input: one subdirectory per dataset key (`README.md`, `chemicals.csv`, …) |
| `workspace_data` | Output: per-key SMILES CSVs, summaries, description markdown |
| `workspace_summary` | Output: aggregated tables (e.g. dataset-of-datasets CSVs) |

Per-dataset pipeline state (detect_files, tags, ess_columns, …) is stored as versioned YAML under `workspace/<key>/` via `VersionedYaml` in `core.py`.

### 2. Secrets — `.env`

Copy [`.env.example`](.env.example) to `.env` next to it (or under `notebooks/` if you prefer) and fill in real values:

```bash
cp .env.example .env
```

Required for LLM steps:

- `OPENAI_API_KEY`
- `OPENAI_API_BASE`
- `OPENAI_MODEL`

Optional Langfuse tracing:

- `LANGFUSE_SECRET_KEY`
- `LANGFUSE_PUBLIC_KEY`
- `LANGFUSE_BASE_URL` (default in example: `http://localhost:3000`)

`core.py` calls `load_dotenv()` on import and initializes the Langfuse client.

## How to run

1. Prepare `workspace/` with one folder per dataset key.
2. Set paths in `paths.yaml`.
3. Configure `.env`.
4. Open notebooks from the `notebooks/` directory (so `from core import *` and local module imports resolve).
5. Run steps **in order** (01 → 05). Early steps write YAML state that later steps read.

Typical flow for a single dataset key (e.g. `ASPECT`):

```text
step01 → detect_files YAML
step02 → tags YAML
step03 → ess_columns YAML (SMILES column mapping)
step04 → workspace_DATA/<key>/{full,unique,broken}_smiles*.csv
step05 → workspace_SUMMARY aggregated tables + templates
```

Most notebooks include a **single-key** cell (set `KEY = "..."`) and an **all-keys** batch cell. Use the single-key cell to debug, then the batch cell for the full workspace.

## Notebooks

| Notebook | Purpose |
|----------|---------|
| **step01_detect_files.ipynb** | Lists files in each dataset folder and uses an LLM to classify them into `data_files`, `description_files`, and `other_files`. Writes `detect_files` YAML per key. |
| **step02_detect_tags.ipynb** | Loads description text, summarizes it, and assigns semantic tags (phase, chemistry domain, regulation, life sciences, materials, etc.). Writes `tags` YAML per key. |
| **step03_detect_SMILES.ipynb** | Inspects candidate CSV files and detects which columns contain SMILES. Writes essential-column / SMILES mapping state (`ess_columns`) used by later processing. |
| **step04_process_SMILES.ipynb** | Reads CSVs via the detected SMILES columns, canonicalizes / splits structures (RDKit), and writes per-key products under `workspace_DATA` (`full_smiles2.csv`, `unique_smiles2.csv`, `broken_smiles2.csv`, `dataset_summary2.json`). Supports parallel batch over all keys. |
| **step05_process_SMILES.ipynb** | Flattens tags and joins unique SMILES with tag templates; builds aggregated summary CSVs under `workspace_SUMMARY` (e.g. `dataset-of-datasets-unique-smiles.csv`, `dataset-of-datasets-templates.csv`). |

## Helper modules

| Module | Role |
|--------|------|
| `core.py` | Env/LLM/Langfuse setup, path loading, `VersionedYaml`, shared config (`LF_CONFIG`, `LLM_BASIC`) |
| `processing.py` | SMILES splitting / `BasicAnalysis` CSV processing |
| `analysis.py` / `structure_utils.py` | RDKit descriptors and structure heuristics |
| `utils.py` | Safe filenames and small utilities |
| `row.py` / `llm_check.py` | Row helpers and LLM connectivity check |

## Notes

- Steps 01–03 call an external LLM; they need network access and valid API credentials.
- Steps 04–05 are largely local RDKit/pandas work once `ess_columns` (and for step 05, `tags`) exist.
- Notebook cell **outputs** may still show paths from earlier runs; path configuration for new runs comes only from `paths.yaml` via `core.py`.
