"""Filter samples from a PLINK dataset via --keep or --remove."""

from __future__ import annotations

from pathlib import Path

from pydantic import BaseModel, Field, model_validator

from statgen_mcp.config import get_config
from statgen_mcp.executors.base import get_executor
from statgen_mcp.utils.file_io import list_output_files, make_run_dir, parse_plink2_log
from statgen_mcp.utils.validation import ValidationError, require_binary, require_file


class FilterSamplesInput(BaseModel):
    input_prefix: str = Field(description="PLINK prefix (without .bed)")
    output_prefix: str = Field(description="Output PLINK prefix")
    keep_file: str | None = Field(
        default=None,
        description=(
            "Path to file of samples to keep (FID IID, whitespace-delimited, no header). "
            "Use for ancestry subsetting or selecting a specific cohort."
        ),
    )
    remove_file: str | None = Field(
        default=None,
        description=(
            "Path to file of samples to remove (FID IID, whitespace-delimited, no header). "
            "Accepts PLINK2 .king.cutoff.out.id directly for kinship removal, "
            "or any FID/IID list (PCA outliers, user-specified, ancestry exclusions)."
        ),
    )
    dry_run: bool = Field(default=False)

    @model_validator(mode="after")
    def at_least_one_file(self) -> "FilterSamplesInput":
        if self.keep_file is None and self.remove_file is None:
            raise ValueError("at least one of keep_file or remove_file is required")
        return self


class FilterSamplesResult(BaseModel):
    status: str
    output_prefix: str = ""
    output_files: list[str] = []
    samples_before: int | None = None
    samples_after: int | None = None
    samples_removed: int | None = None
    log_summary: dict = {}
    error_message: str | None = None
    command_preview: str | None = None
    runtime_seconds: float | None = None


def run_filter_samples(params: FilterSamplesInput) -> FilterSamplesResult:
    cfg = get_config()

    try:
        plink2 = require_binary(cfg.tools.plink2)
        require_file(
            params.input_prefix + ".bed",
            params.input_prefix + ".bim",
            params.input_prefix + ".fam",
        )
        if params.keep_file:
            require_file(params.keep_file)
        if params.remove_file:
            require_file(params.remove_file)
    except ValidationError as e:
        return FilterSamplesResult(status="failed", error_message=str(e))

    run_dir = make_run_dir(cfg.output_dir, "filter_samples")
    out = run_dir / Path(params.output_prefix).name

    cmd = [plink2, "--bfile", params.input_prefix]
    if params.keep_file:
        cmd += ["--keep", params.keep_file]
    if params.remove_file:
        cmd += ["--remove", params.remove_file]
    cmd += ["--make-bed", "--out", str(out)]

    executor = get_executor(cfg.executor.backend)
    job = executor.run(cmd=cmd, job_name="filter_samples", work_dir=run_dir, dry_run=params.dry_run)

    if params.dry_run:
        return FilterSamplesResult(status="dry_run", command_preview=job.stdout)

    if not job.ok:
        return FilterSamplesResult(status="failed", error_message=job.error_message or job.stderr)

    log_path = Path(str(out) + ".log")
    log_summary = parse_plink2_log(log_path)
    output_files = list_output_files(run_dir, [".bed", ".bim", ".fam", ".log"])

    samples_before = log_summary.get("samples")
    samples_after = log_summary.get("samples_after")
    samples_removed = (
        (samples_before - samples_after)
        if samples_before is not None and samples_after is not None
        else None
    )

    return FilterSamplesResult(
        status="success",
        output_prefix=str(out),
        output_files=output_files,
        samples_before=samples_before,
        samples_after=samples_after,
        samples_removed=samples_removed,
        log_summary=log_summary,
        runtime_seconds=job.runtime_seconds,
    )
