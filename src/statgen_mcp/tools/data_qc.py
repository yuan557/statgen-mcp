"""
Data QC tool: MAF, MIND, GENO, HWE filtering via PLINK2.
Accepts VCF, PLINK1 (bed/bim/fam), or PLINK2 native (pgen/psam/pvar) input.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Literal

from pydantic import BaseModel, Field, field_validator, model_validator

from statgen_mcp.config import get_config
from statgen_mcp.executors.base import JobResult, get_executor
from statgen_mcp.utils.file_io import list_output_files, make_run_dir, parse_plink2_log
from statgen_mcp.utils.validation import (
    ValidationError,
    detect_input_format,
    require_binary,
    require_file,
    require_plink_files,
)


class DataQcInput(BaseModel):
    input_prefix: str = Field(
        description="PLINK prefix (without .bed/.pgen extension) or path to VCF/VCF.gz"
    )
    output_prefix: str = Field(
        description="Output PLINK prefix (without extension)"
    )
    maf: float = Field(default=0.01, description="Minor allele frequency threshold (0-0.5)")
    mind: float = Field(default=0.05, description="Per-sample missingness threshold (0-1)")
    geno: float = Field(default=0.05, description="Per-variant missingness threshold (0-1)")
    hwe: float = Field(default=1e-6, description="Hardy-Weinberg p-value threshold")
    hwe_midp: bool = Field(default=True, description="Use mid-p adjustment for HWE test")
    hwe_ctrl_only: bool = Field(
        default=False,
        description=(
            "Apply HWE filter to controls only. "
            "Two defensible approaches: (1) lenient threshold (e.g. 1e-6) with ctrl-only — "
            "protects true disease loci in cases that may deviate from HWE; "
            "(2) stringent threshold (e.g. 5e-8) on the full sample — at that level only "
            "genotyping errors fail, so case inflation is not a practical concern. "
            "Use ctrl-only for standard GWAS QC; full-sample is acceptable if threshold is very stringent. "
            "Ignored when no case-control phenotype is present."
        ),
    )
    biallelic_only: bool = Field(
        default=False,
        description="Keep only biallelic variants (--max-alleles 2). Recommended for KING and PCA.",
    )
    snps_only: bool = Field(
        default=False,
        description="Keep only SNPs, exclude indels (--snps-only just-acgt). Combine with biallelic_only for strictest filtering.",
    )
    variant_source: Literal["wgs", "genotyped", "imputed"] | None = Field(
        default=None,
        description=(
            "Variant source type. 'wgs': no INFO filtering. "
            "'genotyped': keep only directly genotyped variants (array typed flag). "
            "'imputed': keep only variants with INFO score >= min_info_score. "
            "Only applies to VCF input; ignored for PLINK input."
        ),
    )
    min_info_score: float | None = Field(
        default=None,
        description=(
            "Minimum imputation INFO score (required when variant_source='imputed'). "
            "Use >= 0.99 for PCA/kinship input. Use 0.3–0.8 for GWAS input."
        ),
    )
    info_field: str = Field(
        default="R2",
        description="INFO field name for imputation quality (e.g. R2 for Minimac, INFO for IMPUTE2, DR2 for Beagle)",
    )
    typed_flag: str = Field(
        default="TYPED",
        description="INFO flag name for directly genotyped variants (e.g. TYPED for Minimac, GENOTYPED for IMPUTE4)",
    )
    dry_run: bool = Field(default=False, description="Preview command without running")

    @model_validator(mode="after")
    def check_imputed_has_score(self) -> "DataQcInput":
        if self.variant_source == "imputed" and self.min_info_score is None:
            raise ValueError("min_info_score is required when variant_source='imputed'")
        return self

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



class DataQcResult(BaseModel):
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
    variant_source: str | None = None


def _bcftools_filter_expr(params: DataQcInput) -> str | None:
    if params.variant_source == "genotyped":
        return f"INFO/{params.typed_flag}=1"
    if params.variant_source == "imputed":
        return f"INFO/{params.info_field} >= {params.min_info_score}"
    return None


def run_data_qc(params: DataQcInput) -> DataQcResult:
    cfg = get_config()

    try:
        plink2 = cfg.tools.plink2 if params.dry_run else require_binary(cfg.tools.plink2)
    except ValidationError as e:
        return DataQcResult(status="failed", error_message=str(e))

    try:
        fmt = detect_input_format(params.input_prefix)
        if fmt in ("plink1", "plink2"):
            require_plink_files(params.input_prefix, fmt)
        else:
            require_file(params.input_prefix)
    except ValidationError as e:
        return DataQcResult(status="failed", error_message=str(e))

    run_dir = make_run_dir(cfg.output_dir, "data_qc")
    out = run_dir / Path(params.output_prefix).name
    executor = get_executor(cfg.executor.backend)

    warnings = []
    filter_expr = _bcftools_filter_expr(params) if fmt == "vcf" else None

    if params.variant_source and fmt in ("plink1", "plink2"):
        warnings.append("variant_source filtering requires VCF input and was skipped for PLINK input")

    vcf_input = params.input_prefix

    if filter_expr:
        try:
            bcftools = require_binary("bcftools")
        except ValidationError as e:
            return DataQcResult(status="failed", error_message=str(e))

        filtered_vcf = str(run_dir / "filtered.vcf.gz")
        bcftools_cmd = [bcftools, "view", "-i", filter_expr, "-Oz", "-o", filtered_vcf, params.input_prefix]

        if params.dry_run:
            hwe_modifiers = []
            if params.hwe_midp:
                hwe_modifiers.append("midp")
            if params.hwe_ctrl_only:
                hwe_modifiers.append("ctrl-only")
            extra = " ".join(hwe_modifiers)
            plink2_preview = (
                f"plink2 --vcf {filtered_vcf} --maf {params.maf} --mind {params.mind} "
                f"--geno {params.geno} --hwe {params.hwe}{' ' + extra if extra else ''}"
                + (" --max-alleles 2" if params.biallelic_only else "")
                + (" --snps-only just-acgt" if params.snps_only else "")
                + f" --make-bed --out {out}"
            )
            return DataQcResult(
                status="dry_run",
                command_preview=(
                    f"Step 1 (bcftools filter):\n{' '.join(bcftools_cmd)}\n\n"
                    f"Step 2 (PLINK2 QC):\n{plink2_preview}"
                ),
                filters_applied={
                    "variant_source": params.variant_source,
                    "filter_expr": filter_expr,
                    "maf": params.maf,
                    "mind": params.mind,
                    "geno": params.geno,
                    "hwe": params.hwe,
                },
                variant_source=params.variant_source,
            )

        bcftools_job = executor.run(bcftools_cmd, "bcftools_filter", run_dir)
        if not bcftools_job.ok:
            return DataQcResult(status="failed", error_message=bcftools_job.error_message or bcftools_job.stderr)

        vcf_input = filtered_vcf

    cmd = [plink2]
    if fmt == "plink2":
        cmd += ["--pfile", params.input_prefix]
    elif fmt == "plink1":
        cmd += ["--bfile", params.input_prefix]
    else:
        cmd += ["--vcf", vcf_input]

    cmd += [
        "--maf", str(params.maf),
        "--mind", str(params.mind),
        "--geno", str(params.geno),
        "--hwe", str(params.hwe),
    ]
    if params.hwe_midp:
        cmd += ["midp"]
    if params.hwe_ctrl_only:
        cmd += ["ctrl-only"]
    if params.biallelic_only:
        cmd += ["--max-alleles", "2"]
    if params.snps_only:
        cmd += ["--snps-only", "just-acgt"]
    make_flag = "--make-pgen" if fmt == "plink2" else "--make-bed"
    cmd += [make_flag, "--out", str(out)]

    if params.dry_run:
        job: JobResult = executor.run(cmd=cmd, job_name="data_qc", work_dir=run_dir, dry_run=True)
        return DataQcResult(
            status="dry_run",
            command_preview=job.stdout,
            filters_applied={
                "maf": params.maf,
                "mind": params.mind,
                "geno": params.geno,
                "hwe": params.hwe,
            },
            warnings=warnings,
            variant_source=params.variant_source,
        )

    job = executor.run(cmd=cmd, job_name="data_qc", work_dir=run_dir)

    if not job.ok:
        return DataQcResult(status="failed", error_message=job.error_message or job.stderr)

    log_path = Path(str(out) + ".log")
    log_summary = parse_plink2_log(log_path)
    out_exts = [".pgen", ".psam", ".pvar", ".log"] if fmt == "plink2" else [".bed", ".bim", ".fam", ".log"]
    output_files = list_output_files(run_dir, out_exts)

    if log_summary.get("samples_removed_missing", 0) > 0:
        warnings.append(
            f"{log_summary['samples_removed_missing']} samples removed due to high missingness (--mind {params.mind})"
        )

    return DataQcResult(
        status="success",
        output_prefix=str(out),
        output_files=output_files,
        variants_before=log_summary.get("variants_loaded"),
        variants_after=log_summary.get("variants_remaining"),
        samples_before=log_summary.get("samples"),
        log_summary=log_summary,
        filters_applied={
            "variant_source": params.variant_source,
            "maf": params.maf,
            "mind": params.mind,
            "geno": params.geno,
            "hwe": params.hwe,
        },
        warnings=warnings,
        runtime_seconds=job.runtime_seconds,
        variant_source=params.variant_source,
    )
