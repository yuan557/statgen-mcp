"""Input validation helpers shared across tools."""

from __future__ import annotations

import shutil
from pathlib import Path


class ValidationError(Exception):
    pass


def require_file(*paths: str | Path) -> None:
    """Raise ValidationError if any path does not exist."""
    for p in paths:
        if not Path(p).exists():
            raise ValidationError(f"File not found: {p}")


def require_binary(name: str) -> str:
    """Raise ValidationError if binary is not on PATH. Returns resolved path."""
    resolved = shutil.which(name)
    if resolved is None:
        raise ValidationError(
            f"Required binary '{name}' not found on PATH. "
            f"Please install it or set its path in ~/.statgen-mcp/config.toml under [tools]."
        )
    return resolved


def require_range(value: float, lo: float, hi: float, name: str) -> None:
    if not (lo <= value <= hi):
        raise ValidationError(f"{name} must be between {lo} and {hi}, got {value}")


def detect_input_format(prefix_or_vcf: str) -> str:
    """Returns 'plink' or 'vcf' based on file existence."""
    p = Path(prefix_or_vcf)
    if p.suffix in (".vcf", ".gz"):
        return "vcf"
    if Path(str(p) + ".bed").exists():
        return "plink"
    raise ValidationError(
        f"Cannot detect format for '{prefix_or_vcf}'. "
        "Provide a VCF/VCF.gz path or a PLINK prefix (without .bed extension)."
    )
