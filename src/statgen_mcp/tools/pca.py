"""PCA via PLINK2 --pca. Always LD-prunes before running PCA."""

from __future__ import annotations

from pathlib import Path

from pydantic import BaseModel, Field

from statgen_mcp.config import get_config
from statgen_mcp.executors.base import get_executor
from statgen_mcp.utils.file_io import list_output_files, make_run_dir, parse_plink2_log
from statgen_mcp.utils.validation import ValidationError, require_binary, require_file


class PcaInput(BaseModel):
    input_prefix: str = Field(description="PLINK prefix (without .bed)")
    output_prefix: str = Field(description="Output prefix for PCA files")
    n_pcs: int = Field(default=20, description="Number of principal components to compute")
    ld_window: int = Field(default=200, description="LD pruning window size in variants")
    ld_step: int = Field(default=50, description="LD pruning step size")
    ld_r2: float = Field(default=0.2, description="LD pruning r2 threshold")
    approx: bool = Field(
        default=False,
        description="Use --pca approx (faster for very large datasets, slight accuracy tradeoff)",
    )
    high_ld_file: str | None = Field(
        default=None,
        description=(
            "Path to high-LD regions file to exclude before LD pruning (PLINK2 range format: CHR START END NAME). "
            "Strongly recommended for PCA — MHC (chr6:25-35Mb) and inversion regions can drive spurious PCs. "
            "GRCh38 region files are available from Price et al. 2008 / Anderson et al. 2010 or the plinkQC R package."
        ),
    )
    dry_run: bool = Field(default=False)


class PcaResult(BaseModel):
    status: str
    output_prefix: str = ""
    output_files: list[str] = []
    n_pcs: int | None = None
    n_samples: int | None = None
    ld_pruned: bool = True
    log_summary: dict = {}
    error_message: str | None = None
    command_preview: str | None = None
    runtime_seconds: float | None = None


def run_pca(params: PcaInput) -> PcaResult:
    cfg = get_config()

    try:
        plink2 = require_binary(cfg.tools.plink2)
        require_file(
            params.input_prefix + ".bed",
            params.input_prefix + ".bim",
            params.input_prefix + ".fam",
        )
    except ValidationError as e:
        return PcaResult(status="failed", error_message=str(e))

    run_dir = make_run_dir(cfg.output_dir, "pca")
    out = run_dir / Path(params.output_prefix).name
    executor = get_executor(cfg.executor.backend)
    pruned_prefix = str(out) + "_pruned"

    prune_cmd = [
        plink2,
        "--bfile", params.input_prefix,
        "--indep-pairwise",
        str(params.ld_window),
        str(params.ld_step),
        str(params.ld_r2),
    ]
    if params.high_ld_file:
        prune_cmd += ["--exclude", "range", params.high_ld_file]
    prune_cmd += ["--out", str(out) + "_prune"]

    if params.dry_run:
        ld_note = f" (excluding high-LD regions: {params.high_ld_file})" if params.high_ld_file else ""
        return PcaResult(
            status="dry_run",
            command_preview=(
                f"Step 1 (LD prune{ld_note}):\n{' '.join(str(c) for c in prune_cmd)}\n\n"
                f"Step 2 (PCA): plink2 --bfile {pruned_prefix} --pca {params.n_pcs} --out {out}"
            ),
            n_pcs=params.n_pcs,
        )

    prune_job = executor.run(prune_cmd, "ld_prune", run_dir)
    if not prune_job.ok:
        return PcaResult(status="failed", error_message=prune_job.error_message or prune_job.stderr)

    extract_cmd = [
        plink2,
        "--bfile", params.input_prefix,
        "--extract", str(out) + "_prune.prune.in",
        "--make-bed",
        "--out", pruned_prefix,
    ]
    extract_job = executor.run(extract_cmd, "ld_extract", run_dir)
    if not extract_job.ok:
        return PcaResult(status="failed", error_message=extract_job.error_message or extract_job.stderr)

    pca_input = pruned_prefix

    pca_cmd = [
        plink2,
        "--bfile", pca_input,
        "--pca", str(params.n_pcs),
    ]
    if params.approx:
        pca_cmd.append("approx")
    pca_cmd += ["--out", str(out)]

    pca_job = executor.run(pca_cmd, "pca", run_dir)
    if not pca_job.ok:
        return PcaResult(status="failed", error_message=pca_job.error_message or pca_job.stderr)

    log_path = Path(str(out) + ".log")
    log_summary = parse_plink2_log(log_path)
    output_files = list_output_files(run_dir, [".eigenvec", ".eigenval", ".log"])

    return PcaResult(
        status="success",
        output_prefix=str(out),
        output_files=output_files,
        n_pcs=params.n_pcs,
        n_samples=log_summary.get("samples"),
        ld_pruned=True,
        log_summary=log_summary,
        runtime_seconds=pca_job.runtime_seconds,
    )
