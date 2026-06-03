"""
Relatedness filtering via PLINK2 --king-cutoff.
Removes one sample from each pair exceeding the kinship threshold.
"""

from __future__ import annotations

from pathlib import Path

from pydantic import BaseModel, Field

from statgen_mcp.config import get_config
from statgen_mcp.executors.base import JobResult, get_executor
from statgen_mcp.utils.file_io import list_output_files, make_run_dir, parse_plink2_log
from statgen_mcp.utils.validation import ValidationError, require_binary, require_file


KINSHIP_THRESHOLDS = {
    "first_degree": 0.177,
    "second_degree": 0.0884,
    "third_degree": 0.0442,
    "unrelated": 0.0221,
}


class KingInput(BaseModel):
    input_prefix: str = Field(description="PLINK prefix (without .bed)")
    output_prefix: str = Field(description="Output PLINK prefix")
    king_cutoff: float = Field(
        default=0.0442,
        description=(
            "Kinship coefficient cutoff. One sample from each pair above this is removed. "
            "Common values: 0.177 (1st degree), 0.0884 (2nd), 0.0442 (3rd / typical GWAS cutoff)"
        ),
    )
    dry_run: bool = Field(default=False)


class KingResult(BaseModel):
    status: str
    output_prefix: str = ""
    output_files: list[str] = []
    samples_before: int | None = None
    kinship_cutoff_used: float | None = None
    log_summary: dict = {}
    error_message: str | None = None
    command_preview: str | None = None
    runtime_seconds: float | None = None


def run_king_cutoff(params: KingInput) -> KingResult:
    cfg = get_config()

    try:
        plink2 = require_binary(cfg.tools.plink2)
        require_file(
            params.input_prefix + ".bed",
            params.input_prefix + ".bim",
            params.input_prefix + ".fam",
        )
    except ValidationError as e:
        return KingResult(status="failed", error_message=str(e))

    run_dir = make_run_dir(cfg.output_dir, "king_cutoff")
    out = run_dir / Path(params.output_prefix).name

    cmd = [
        plink2,
        "--bfile", params.input_prefix,
        "--king-cutoff", str(params.king_cutoff),
        "--make-bed",
        "--out", str(out),
    ]

    executor = get_executor(cfg.executor.backend)
    job: JobResult = executor.run(
        cmd=cmd,
        job_name="king_cutoff",
        work_dir=run_dir,
        dry_run=params.dry_run,
    )

    if params.dry_run:
        return KingResult(
            status="dry_run",
            command_preview=job.stdout,
            kinship_cutoff_used=params.king_cutoff,
        )

    if not job.ok:
        return KingResult(status="failed", error_message=job.error_message or job.stderr)

    log_path = Path(str(out) + ".log")
    log_summary = parse_plink2_log(log_path)
    output_files = list_output_files(
        run_dir,
        [".bed", ".bim", ".fam", ".log", ".king.cutoff.out.id", ".king.cutoff.in.id"]
    )

    return KingResult(
        status="success",
        output_prefix=str(out),
        output_files=output_files,
        samples_before=log_summary.get("samples"),
        log_summary=log_summary,
        kinship_cutoff_used=params.king_cutoff,
        runtime_seconds=job.runtime_seconds,
    )
