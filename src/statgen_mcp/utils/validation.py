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
    """Returns 'vcf', 'plink2' (pgen/psam/pvar), or 'plink1' (bed/bim/fam)."""
    p = Path(prefix_or_vcf)
    if p.suffix in (".vcf", ".gz"):
        return "vcf"
    if Path(str(p) + ".pgen").exists():
        return "plink2"
    if Path(str(p) + ".bed").exists():
        return "plink1"
    raise ValidationError(
        f"Cannot detect format for '{prefix_or_vcf}'. "
        "Provide a VCF/VCF.gz path or a PLINK prefix (without .bed/.pgen extension)."
    )


def require_plink_files(prefix: str, fmt: str) -> None:
    """Validate that all expected PLINK binary files exist for the given format."""
    if fmt == "plink2":
        require_file(prefix + ".pgen", prefix + ".psam", prefix + ".pvar")
    else:
        require_file(prefix + ".bed", prefix + ".bim", prefix + ".fam")
