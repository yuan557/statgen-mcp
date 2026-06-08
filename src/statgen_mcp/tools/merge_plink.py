"""Merge multiple per-chromosome PLINK datasets into one genome-wide file."""

from __future__ import annotations

from pathlib import Path

from pydantic import BaseModel, Field, field_validator

from statgen_mcp.config import get_config
from statgen_mcp.executors.base import get_executor
from statgen_mcp.utils.file_io import list_output_files, make_run_dir, parse_plink2_log
from statgen_mcp.utils.validation import (
    ValidationError,
    detect_input_format,
    require_binary,
    require_plink_files,
)


class MergePlinkInput(BaseModel):
    input_prefixes: list[str] = Field(
        description="List of PLINK prefixes to merge (e.g. one per chromosome)"
    )
    output_prefix: str = Field(
        description="Output PLINK prefix (without extension)"
    )
    slurm_opts: dict = Field(default_factory=dict, description="Override SLURM parameters for this job")
    dry_run: bool = Field(default=False, description="Preview command without running")

    @field_validator("input_prefixes")
    @classmethod
    def check_at_least_two(cls, v: list[str]) -> list[str]:
        if len(v) < 2:
            raise ValueError("input_prefixes must contain at least two prefixes to merge")
        return v


class MergePlinkResult(BaseModel):
    status: str
    output_prefix: str = ""
    output_files: list[str] = []
    n_files_merged: int = 0
    log_summary: dict = {}
    error_message: str | None = None
    command_preview: str | None = None
    runtime_seconds: float | None = None


def run_merge_plink(params: MergePlinkInput) -> MergePlinkResult:
    cfg = get_config()

    try:
        plink2 = cfg.tools.plink2 if params.dry_run else require_binary(cfg.tools.plink2)
    except ValidationError as e:
        return MergePlinkResult(status="failed", error_message=str(e))

    try:
        fmt = detect_input_format(params.input_prefixes[0])
        if fmt not in ("plink1", "plink2"):
            raise ValidationError("merge_plink requires PLINK input (bed/bim/fam or pgen/psam/pvar), not VCF")
        for prefix in params.input_prefixes:
            require_plink_files(prefix, fmt)
    except ValidationError as e:
        return MergePlinkResult(status="failed", error_message=str(e))

    run_dir = make_run_dir(cfg.output_dir, "merge_plink")
    out = run_dir / Path(params.output_prefix).name
    executor = get_executor(cfg.executor.backend)

    merge_list = run_dir / "merge_list.txt"
    merge_list.write_text("\n".join(params.input_prefixes) + "\n")

    fmt_flag = "pfile" if fmt == "plink2" else "bfile"
    make_flag = "--make-pgen" if fmt == "plink2" else "--make-bed"

    cmd = [
        plink2,
        "--pmerge-list", str(merge_list), fmt_flag,
        make_flag, "--out", str(out),
    ]

    if params.dry_run:
        job = executor.run(cmd=cmd, job_name="merge_plink", work_dir=run_dir, dry_run=True)
        return MergePlinkResult(
            status="dry_run",
            command_preview=job.stdout,
            n_files_merged=len(params.input_prefixes),
        )

    job = executor.run(cmd=cmd, job_name="merge_plink", work_dir=run_dir, slurm_opts=params.slurm_opts)

    if not job.ok:
        return MergePlinkResult(status="failed", error_message=job.error_message or job.stderr)

    log_path = Path(str(out) + ".log")
    log_summary = parse_plink2_log(log_path)
    out_exts = [".pgen", ".psam", ".pvar", ".log"] if fmt == "plink2" else [".bed", ".bim", ".fam", ".log"]
    output_files = list_output_files(run_dir, out_exts)

    return MergePlinkResult(
        status="success",
        output_prefix=str(out),
        output_files=output_files,
        n_files_merged=len(params.input_prefixes),
        log_summary=log_summary,
        runtime_seconds=job.runtime_seconds,
    )
