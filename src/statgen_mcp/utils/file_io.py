"""File helpers: run directory creation, log parsing, PLINK output parsing."""

from __future__ import annotations

import re
from datetime import datetime
from pathlib import Path


def make_run_dir(base: str, job_name: str) -> Path:
    """Create a timestamped run directory and return it."""
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = Path(base).expanduser() / f"{job_name}_{ts}"
    run_dir.mkdir(parents=True, exist_ok=True)
    return run_dir


def list_output_files(directory: Path, extensions: list[str]) -> list[str]:
    """Return files in directory matching given extensions."""
    found = []
    for ext in extensions:
        found.extend(str(f) for f in sorted(directory.glob(f"*{ext}")))
    return found


def parse_plink2_log(log_path: Path) -> dict:
    """Parse a PLINK2 .log file and extract key numbers."""
    if not log_path.exists():
        return {}

    text = log_path.read_text()
    result: dict = {}

    m = re.search(r"(\d+) sample[s]? \(", text)
    if m:
        result["samples"] = int(m.group(1))

    m = re.search(r"(\d+) variant[s]? loaded", text)
    if m:
        result["variants_loaded"] = int(m.group(1))

    m = re.search(r"(\d+) variant[s]? remaining", text)
    if m:
        result["variants_remaining"] = int(m.group(1))

    m = re.search(r"(\d+) sample[s]?.*remaining after", text)
    if m:
        result["samples_after"] = int(m.group(1))

    m = re.search(r"(\d+) sample[s]? removed due to missing", text)
    if m:
        result["samples_removed_missing"] = int(m.group(1))

    m = re.search(r"(\d+) variant[s]? removed due to Hardy-Weinberg", text)
    if m:
        result["variants_removed_hwe"] = int(m.group(1))

    m = re.search(r"(\d+) variant[s]? removed due to allele frequency", text)
    if m:
        result["variants_removed_maf"] = int(m.group(1))

    if re.search(r"linear regression", text, re.IGNORECASE):
        result["regression_type"] = "quantitative"
    elif re.search(r"logistic", text, re.IGNORECASE):
        result["regression_type"] = "binary"

    return result
