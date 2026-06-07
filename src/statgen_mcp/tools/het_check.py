"""Heterozygosity QC via PLINK2 --het. Flags contaminated and inbred samples."""

from __future__ import annotations

import csv
from pathlib import Path

from pydantic import BaseModel, Field

from statgen_mcp.config import get_config
from statgen_mcp.executors.base import get_executor
from statgen_mcp.utils.file_io import make_run_dir
from statgen_mcp.utils.validation import (
    ValidationError,
    detect_input_format,
    require_binary,
    require_plink_files,
)


class HetCheckInput(BaseModel):
    input_prefix: str = Field(description="PLINK prefix (without .bed/.pgen extension)")
    output_prefix: str = Field(description="Output prefix for het file and outlier list")
    sd_threshold: float = Field(
        default=3.0,
        description=(
            "Flag samples whose F statistic deviates more than this many SD from the mean. "
            "Standard is 3 SD. F > mean+threshold = deficient heterozygosity (inbreeding). "
            "F < mean-threshold = excess heterozygosity (contamination)."
        ),
    )
    write_remove_file: bool = Field(
        default=True,
        description="Write a FID/IID remove file of outliers, ready for filter_samples.",
    )
    dry_run: bool = Field(default=False)


class HetOutlier(BaseModel):
    fid: str
    iid: str
    f_stat: float
    direction: str


class HetCheckResult(BaseModel):
    status: str
    output_prefix: str = ""
    n_samples: int | None = None
    n_outliers: int | None = None
    mean_f: float | None = None
    sd_f: float | None = None
    outliers: list[HetOutlier] = []
    remove_file: str | None = None
    het_file: str | None = None
    error_message: str | None = None
    command_preview: str | None = None
    runtime_seconds: float | None = None


def _parse_het_file(path: Path, sd_threshold: float) -> tuple[list[HetOutlier], float, float, int]:
    import math

    rows = []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            fid = row.get("FID") or row.get("#FID", "")
            iid = row.get("IID", "")
            try:
                f_stat = float(row["F"])
                rows.append((fid, iid, f_stat))
            except (KeyError, ValueError):
                continue

    if not rows:
        return [], 0.0, 0.0, 0

    f_vals = [r[2] for r in rows]
    n = len(f_vals)
    mean_f = sum(f_vals) / n
    sd_f = math.sqrt(sum((x - mean_f) ** 2 for x in f_vals) / n)

    outliers = []
    for fid, iid, f_stat in rows:
        if sd_f == 0:
            continue
        z = (f_stat - mean_f) / sd_f
        if z > sd_threshold:
            outliers.append(HetOutlier(fid=fid, iid=iid, f_stat=round(f_stat, 6), direction="low_het"))
        elif z < -sd_threshold:
            outliers.append(HetOutlier(fid=fid, iid=iid, f_stat=round(f_stat, 6), direction="excess_het"))

    return outliers, round(mean_f, 6), round(sd_f, 6), n


def run_het_check(params: HetCheckInput) -> HetCheckResult:
    cfg = get_config()

    try:
        plink2 = cfg.tools.plink2 if params.dry_run else require_binary(cfg.tools.plink2)
        fmt = detect_input_format(params.input_prefix)
        require_plink_files(params.input_prefix, fmt)
    except ValidationError as e:
        return HetCheckResult(status="failed", error_message=str(e))

    run_dir = make_run_dir(cfg.output_dir, "het_check")
    out = run_dir / Path(params.output_prefix).name

    plink_flag = "--pfile" if fmt == "plink2" else "--bfile"
    cmd = [plink2, plink_flag, params.input_prefix, "--het", "--out", str(out)]

    executor = get_executor(cfg.executor.backend)
    job = executor.run(cmd=cmd, job_name="het_check", work_dir=run_dir, dry_run=params.dry_run)

    if params.dry_run:
        return HetCheckResult(status="dry_run", command_preview=job.stdout)

    if not job.ok:
        return HetCheckResult(status="failed", error_message=job.error_message or job.stderr)

    het_path = Path(str(out) + ".het")
    if not het_path.exists():
        return HetCheckResult(status="failed", error_message=f"Expected het file not found: {het_path}")

    outliers, mean_f, sd_f, n_samples = _parse_het_file(het_path, params.sd_threshold)

    remove_file = None
    if params.write_remove_file and outliers:
        remove_path = run_dir / "het_outliers_remove.txt"
        with open(remove_path, "w") as f:
            for o in outliers:
                f.write(f"{o.fid}\t{o.iid}\n")
        remove_file = str(remove_path)

    return HetCheckResult(
        status="success",
        output_prefix=str(out),
        n_samples=n_samples,
        n_outliers=len(outliers),
        mean_f=mean_f,
        sd_f=sd_f,
        outliers=outliers,
        remove_file=remove_file,
        het_file=str(het_path),
        runtime_seconds=job.runtime_seconds,
    )
