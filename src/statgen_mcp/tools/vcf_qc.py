"""
VCF QC tool: MAF, MIND, GENO, HWE filtering via PLINK2.
Also handles VCF -> PLINK conversion if needed.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from pydantic import BaseModel, Field, field_validator

from statgen_mcp.config import get_config
from statgen_mcp.executors.base import JobResult, get_executor
from statgen_mcp.utils.file_io import list_output_files, make_run_dir, parse_plink2_log
from statgen_mcp.utils.validation import (
    ValidationError,
    detect_input_format,
    require_binary,
    require_file,
)


class VcfQcInput(BaseModel):
    input_prefix: str = Field(
        description="PLINK prefix (without .bed) or path to VCF/VCF.gz"
    )
    output_prefix: str = Field(
        description="Output PLINK prefix (without extension)"
    )
    maf: float = Field(default=0.01, description="Minor allele frequency threshold (0-0.5)")
    mind: float = Field(default=0.05, description="Per-sample missingness threshold (0-1)")
    geno: float = Field(default=0.05, description="Per-variant missingness threshold (0-1)")
    hwe: float = Field(default=1e-6, description="Hardy-Weinberg p-value threshold")
    hwe_midp: bool = Field(default=True, description="Use mid-p adjustment for HWE test")
    dry_run: bool = Field(default=False, description="Preview command without running")

    @field_validator("maf")
    @classmethod
    def check_maf(cls, v: float) -> float:
        if not (0 <= v <= 0.5):
            raise ValueError("maf must be between 0 and 0.5")
        return v

    @field_validator("mind", "geno")
    @classmethod
    def check_missingness(cls, v: float) -> float:
        if not (0 <= v <= 1):
            raise ValueError("must be between 0 and 1")
        return v


class VcfQcResult(BaseModel):
    status: str
    output_prefix: str = ""
    output_files: list[str] = []
    samples_before: int | None = None
    variants_before: int | None = None
    variants_after: int | None = None
    filters_applied: dict[str, Any] = {}
    warnings: list[str] = []
    log_summary: dict = {}
    error_message: str | None = None
    command_preview: str | None = None
    runtime_seconds: float | None = None


def run_vcf_qc(params: VcfQcInput) -> VcfQcResult:
    cfg = get_config()

    try:
        plink2 = require_binary(cfg.tools.plink2)
    except ValidationError as e:
        return VcfQcResult(status="failed", error_message=str(e))

    try:
        fmt = detect_input_format(params.input_prefix)
        if fmt == "plink":
            require_file(
                params.input_prefix + ".bed",
                params.input_prefix + ".bim",
                params.input_prefix + ".fam",
            )
        else:
            require_file(params.input_prefix)
    except ValidationError as e:
        return VcfQcResult(status="failed", error_message=str(e))

    run_dir = make_run_dir(cfg.output_dir, "vcf_qc")
    out = run_dir / Path(params.output_prefix).name

    cmd = [plink2]

    if fmt == "plink":
        cmd += ["--bfile", params.input_prefix]
    else:
        cmd += ["--vcf", params.input_prefix]

    cmd += [
        "--maf", str(params.maf),
        "--mind", str(params.mind),
        "--geno", str(params.geno),
        "--hwe", str(params.hwe),
    ]

    if params.hwe_midp:
        cmd += ["midp"]

    cmd += ["--make-bed", "--out", str(out)]

    executor = get_executor(cfg.executor.backend)
    job: JobResult = executor.run(
        cmd=cmd,
        job_name="vcf_qc",
        work_dir=run_dir,
        dry_run=params.dry_run,
    )

    if params.dry_run:
        return VcfQcResult(
            status="dry_run",
            command_preview=job.stdout,
            filters_applied={
                "maf": params.maf,
                "mind": params.mind,
                "geno": params.geno,
                "hwe": params.hwe,
            },
        )

    if not job.ok:
        return VcfQcResult(
            status="failed",
            error_message=job.error_message or job.stderr,
        )

    log_path = Path(str(out) + ".log")
    log_summary = parse_plink2_log(log_path)
    output_files = list_output_files(run_dir, [".bed", ".bim", ".fam", ".log"])

    warnings = []
    if log_summary.get("samples_removed_missing", 0) > 0:
        warnings.append(
            f"{log_summary['samples_removed_missing']} samples removed due to high missingness (--mind {params.mind})"
        )

    return VcfQcResult(
        status="success",
        output_prefix=str(out),
        output_files=output_files,
        variants_before=log_summary.get("variants_loaded"),
        variants_after=log_summary.get("variants_remaining"),
        samples_before=log_summary.get("samples"),
        log_summary=log_summary,
        filters_applied={
            "maf": params.maf,
            "mind": params.mind,
            "geno": params.geno,
            "hwe": params.hwe,
        },
        warnings=warnings,
        runtime_seconds=job.runtime_seconds,
    )
