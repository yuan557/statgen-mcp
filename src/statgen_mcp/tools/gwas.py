"""GWAS via PLINK2 --glm. Handles covariate assembly from user file + PCA eigenvec."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Literal

from pydantic import BaseModel, Field

from statgen_mcp.config import get_config
from statgen_mcp.executors.base import get_executor
from statgen_mcp.utils.file_io import list_output_files, make_run_dir, parse_plink2_log
from statgen_mcp.utils.validation import (
    ValidationError,
    detect_input_format,
    require_binary,
    require_file,
    require_plink_files,
)


class GwasInput(BaseModel):
    input_prefix: str = Field(description="QC'd PLINK prefix (without .bed/.pgen extension)")
    output_prefix: str = Field(description="Output prefix for GWAS results")
    pheno_file: str = Field(
        description="Phenotype file (FID IID TRAIT1 [TRAIT2 ...]). Space/tab-delimited."
    )
    pheno_name: str | None = Field(
        default=None,
        description="Column name of the trait to test. Required if pheno_file has multiple trait columns.",
    )
    covar_file: str | None = Field(
        default=None,
        description="Covariate file (FID IID COV1 COV2 ...). Merged with PCA eigenvec if both provided.",
    )
    eigenvec_file: str | None = Field(
        default=None,
        description="PCA eigenvec file from run_pca. First n_pcs PCs are added as covariates.",
    )
    n_pcs: int = Field(
        default=10,
        description="Number of PCs to include as covariates (1-20). Claude should choose based on eigenvalues.",
    )
    trait_type: Literal["binary", "quantitative"] | None = Field(
        default=None,
        description="Trait type. None = auto-detect from phenotype values. Binary triggers firth-fallback.",
    )
    firth_fallback: bool = Field(
        default=True,
        description="Use Firth regression as fallback for binary traits when standard logistic fails.",
    )
    ci: bool = Field(default=True, description="Include 95% confidence intervals in output.")
    dry_run: bool = Field(default=False)


class GwasResult(BaseModel):
    status: str
    output_prefix: str = ""
    output_files: list[str] = []
    n_samples: int | None = None
    n_variants_tested: int | None = None
    trait_type_used: str | None = None
    covar_file_used: str | None = None
    log_summary: dict = {}
    error_message: str | None = None
    command_preview: str | None = None
    runtime_seconds: float | None = None


def _build_covar_file(
    covar_file: str | None,
    eigenvec_file: str | None,
    n_pcs: int,
    run_dir: Path,
) -> str | None:
    """Merge user covar file and PCA eigenvec into a single covariate file."""
    if not covar_file and not eigenvec_file:
        return None

    # Load user covariates keyed by (FID, IID)
    covar_rows: dict[tuple[str, str], dict[str, str]] = {}
    covar_cols: list[str] = []
    if covar_file:
        with open(covar_file) as f:
            reader = csv.DictReader(f, delimiter="\t")
            if reader.fieldnames is None:
                # try space-delimited
                f.seek(0)
                reader = csv.DictReader(f, delimiter=" ")
            for row in reader:
                fid = row.get("FID") or row.get("#FID", "")
                iid = row.get("IID", "")
                covar_rows[(fid, iid)] = row
        covar_cols = [
            c for c in (reader.fieldnames or [])
            if c not in ("#FID", "FID", "IID")
        ]

    # Load PCA eigenvec
    pc_cols: list[str] = []
    pc_rows: dict[tuple[str, str], dict[str, str]] = {}
    if eigenvec_file:
        with open(eigenvec_file) as f:
            reader = csv.DictReader(f, delimiter="\t")
            if reader.fieldnames is None:
                f.seek(0)
                reader = csv.DictReader(f, delimiter=" ")
            all_pc_cols = [
                c for c in (reader.fieldnames or [])
                if c not in ("#FID", "FID", "IID")
            ]
            pc_cols = all_pc_cols[:n_pcs]
            for row in reader:
                fid = row.get("FID") or row.get("#FID", "")
                iid = row.get("IID", "")
                pc_rows[(fid, iid)] = row

    # Determine sample universe — intersection of both sources if both provided
    if covar_rows and pc_rows:
        samples = [k for k in covar_rows if k in pc_rows]
    elif covar_rows:
        samples = list(covar_rows.keys())
    else:
        samples = list(pc_rows.keys())

    out_path = run_dir / "combined_covariates.txt"
    all_cols = covar_cols + pc_cols

    with open(out_path, "w") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["#FID", "IID"] + all_cols)
        for fid, iid in samples:
            row_out = [fid, iid]
            for col in covar_cols:
                row_out.append(covar_rows[(fid, iid)].get(col, "NA"))
            for col in pc_cols:
                row_out.append(pc_rows[(fid, iid)].get(col, "NA"))
            writer.writerow(row_out)

    return str(out_path)


def run_gwas(params: GwasInput) -> GwasResult:
    cfg = get_config()

    try:
        plink2 = cfg.tools.plink2 if params.dry_run else require_binary(cfg.tools.plink2)
        fmt = detect_input_format(params.input_prefix)
        require_plink_files(params.input_prefix, fmt)
        require_file(params.pheno_file)
        if params.covar_file:
            require_file(params.covar_file)
        if params.eigenvec_file:
            require_file(params.eigenvec_file)
    except ValidationError as e:
        return GwasResult(status="failed", error_message=str(e))

    run_dir = make_run_dir(cfg.output_dir, "gwas")
    out = run_dir / Path(params.output_prefix).name

    covar_path = _build_covar_file(
        params.covar_file, params.eigenvec_file, params.n_pcs, run_dir
    )

    glm_modifiers = ["hide-covar"]
    if params.firth_fallback and params.trait_type != "quantitative":
        glm_modifiers.append("firth-fallback")

    plink_flag = "--pfile" if fmt == "plink2" else "--bfile"

    cmd = [
        plink2,
        plink_flag, params.input_prefix,
        "--no-psam-pheno",
        "--pheno", params.pheno_file,
    ]
    if params.pheno_name:
        cmd += ["--pheno-name", params.pheno_name]
    if covar_path:
        cmd += ["--covar", covar_path, "--covar-variance-standardize"]
    cmd += ["--glm"] + glm_modifiers
    if params.ci:
        cmd += ["--ci", "0.95"]
    cmd += ["--out", str(out)]

    executor = get_executor(cfg.executor.backend)
    job = executor.run(cmd=cmd, job_name="gwas", work_dir=run_dir, dry_run=params.dry_run)

    if params.dry_run:
        return GwasResult(
            status="dry_run",
            command_preview=job.stdout,
            covar_file_used=covar_path,
        )

    if not job.ok:
        return GwasResult(status="failed", error_message=job.error_message or job.stderr)

    log_path = Path(str(out) + ".log")
    log_summary = parse_plink2_log(log_path)
    output_files = list_output_files(run_dir, [".glm.linear", ".glm.logistic", ".glm.logistic.hybrid", ".log"])

    trait_type_used = _detect_trait_type(log_summary)

    return GwasResult(
        status="success",
        output_prefix=str(out),
        output_files=output_files,
        n_samples=log_summary.get("samples_after") or log_summary.get("samples"),
        n_variants_tested=log_summary.get("variants_remaining") or log_summary.get("variants_loaded"),
        trait_type_used=trait_type_used,
        covar_file_used=covar_path,
        log_summary=log_summary,
        runtime_seconds=job.runtime_seconds,
    )


def _detect_trait_type(log_summary: dict) -> str | None:
    raw = log_summary.get("regression_type")
    if raw:
        return raw
    return None
