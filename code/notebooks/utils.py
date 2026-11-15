import re
import unicodedata

def safe_filename(value: str, default: str = "untitled", max_length: int = 255) -> str:
    """
    Convert an arbitrary string into a safe *nix filename.
    Rules:
    - Normalize unicode to ASCII when possible.
    - Replace path separators and whitespace with underscores.
    - Keep only [A-Za-z0-9._-].
    - Collapse repeated separators.
    - Avoid empty names and reserved dot names ('.', '..').
    - Trim to max_length (default: 255 bytes/chars).
    """
    if max_length < 1:
        raise ValueError("max_length must be >= 1")
    if value is None:
        value = ""
    elif not isinstance(value, str):
        value = str(value)
    # Normalize to a predictable ASCII subset.
    normalized = unicodedata.normalize("NFKD", value)
    normalized = normalized.encode("ascii", "ignore").decode("ascii")
    # Remove control chars and null bytes.
    normalized = "".join(ch for ch in normalized if ch >= " " and ch != "\x7f")
    # Replace risky separators and whitespace.
    normalized = normalized.replace("/", "_").replace("\\", "_")
    normalized = re.sub(r"\s+", "_", normalized)
    # Keep a conservative set of characters.
    normalized = re.sub(r"[^A-Za-z0-9._-]", "_", normalized)
    # Collapse repeated separators.
    normalized = re.sub(r"_+", "_", normalized)
    normalized = re.sub(r"-{2,}", "-", normalized)
    normalized = re.sub(r"\.{2,}", ".", normalized)
    # Strip leading/trailing separators.
    candidate = normalized.strip("._-")
    # Avoid invalid/special dot-only names.
    if candidate in {"", ".", ".."}:
        candidate = default
    # Ensure fallback is safe too.
    candidate = re.sub(r"[^A-Za-z0-9._-]", "_", candidate)
    candidate = candidate.strip("._-") or "untitled"
    if len(candidate) > max_length:
        candidate = candidate[:max_length].rstrip("._-") or "untitled"
    return candidate
