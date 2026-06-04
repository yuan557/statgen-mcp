"""
Executor abstraction.
LocalExecutor runs subprocesses directly.
SlurmExecutor submits sbatch jobs, polls until completion, and returns results.
"""

from __future__ import annotations

import re
import subprocess
import time
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class JobResult:
    status: str                        # "success" | "failed" | "dry_run"
    returncode: int = 0
    stdout: str = ""
    stderr: str = ""
    output_files: list[str] = field(default_factory=list)
    runtime_seconds: float = 0.0
    job_id: str | None = None
    error_message: str | None = None

    @property
    def ok(self) -> bool:
        return self.status == "success"


class BaseExecutor(ABC):
    @abstractmethod
    def run(
        self,
        cmd: list[str],
        job_name: str,
        work_dir: Path,
        cpus: int = 1,
        memory_gb: int = 4,
        dry_run: bool = False,
        slurm_opts: dict | None = None,
    ) -> JobResult:
        ...

    def preview(self, cmd: list[str]) -> str:
        return " ".join(str(c) for c in cmd)


class LocalExecutor(BaseExecutor):
    """Runs commands locally via subprocess."""

    def run(
        self,
        cmd: list[str],
        job_name: str,
        work_dir: Path,
        cpus: int = 1,
        memory_gb: int = 4,
        dry_run: bool = False,
        slurm_opts: dict | None = None,
    ) -> JobResult:
        work_dir.mkdir(parents=True, exist_ok=True)

        if dry_run:
            return JobResult(
                status="dry_run",
                stdout=f"[DRY RUN] Would execute:\n{self.preview(cmd)}",
            )

        start = time.time()
        try:
            result = subprocess.run(
                cmd,
                cwd=work_dir,
                capture_output=True,
                text=True,
            )
            elapsed = time.time() - start

            if result.returncode == 0:
                return JobResult(
                    status="success",
                    returncode=result.returncode,
                    stdout=result.stdout,
                    stderr=result.stderr,
                    runtime_seconds=elapsed,
                )
            else:
                return JobResult(
                    status="failed",
                    returncode=result.returncode,
                    stdout=result.stdout,
                    stderr=result.stderr,
                    runtime_seconds=elapsed,
                    error_message=result.stderr.strip().splitlines()[-1] if result.stderr else "Non-zero exit",
                )

        except FileNotFoundError as e:
            return JobResult(
                status="failed",
                error_message=f"Binary not found: {e}. Is it installed and on PATH?",
            )
        except Exception as e:
            return JobResult(
                status="failed",
                error_message=str(e),
            )


class SlurmExecutor(BaseExecutor):
    """Submits jobs via sbatch, polls until completion, returns results."""

    def run(
        self,
        cmd: list[str],
        job_name: str,
        work_dir: Path,
        cpus: int = 1,
        memory_gb: int = 4,
        dry_run: bool = False,
        slurm_opts: dict | None = None,
    ) -> JobResult:
        from statgen_mcp.config import get_config
        cfg = get_config()
        slurm_cfg = cfg.slurm
        opts = slurm_opts or {}

        work_dir.mkdir(parents=True, exist_ok=True)

        partition  = opts.get("partition",  slurm_cfg.partition)
        time_limit = opts.get("time",       slurm_cfg.time)
        memory     = opts.get("mem",        slurm_cfg.default_memory)
        ncpus      = opts.get("cpus",       max(cpus, slurm_cfg.default_cpus))
        account    = opts.get("account",    slurm_cfg.account)

        stdout_file = work_dir / f"{job_name}.out"
        stderr_file = work_dir / f"{job_name}.err"
        script_file = work_dir / f"{job_name}.sh"

        script_lines = [
            "#!/bin/bash",
            f"#SBATCH --job-name={job_name}",
            f"#SBATCH --output={stdout_file}",
            f"#SBATCH --error={stderr_file}",
            f"#SBATCH --partition={partition}",
            f"#SBATCH --time={time_limit}",
            f"#SBATCH --mem={memory}",
            f"#SBATCH --cpus-per-task={ncpus}",
        ]
        if account:
            script_lines.append(f"#SBATCH --account={account}")
        script_lines += ["", self.preview(cmd), ""]
        script_file.write_text("\n".join(script_lines))

        if dry_run:
            return JobResult(
                status="dry_run",
                stdout=(
                    f"[DRY RUN] Would submit:\n{script_file}\n\n"
                    f"--- script ---\n{''.join(script_lines)}\n\n"
                    f"--- command ---\n{self.preview(cmd)}"
                ),
            )

        start = time.time()

        # Submit
        try:
            submit = subprocess.run(
                ["sbatch", str(script_file)],
                capture_output=True, text=True, cwd=work_dir,
            )
        except FileNotFoundError:
            return JobResult(status="failed", error_message="sbatch not found — is SLURM available?")

        if submit.returncode != 0:
            return JobResult(status="failed", error_message=submit.stderr.strip())

        m = re.search(r"(\d+)", submit.stdout)
        if not m:
            return JobResult(status="failed", error_message=f"Could not parse job ID from: {submit.stdout!r}")
        job_id = m.group(1)

        # Poll until job leaves the queue
        waited = 0
        while waited < slurm_cfg.max_wait:
            time.sleep(slurm_cfg.poll_interval)
            waited += slurm_cfg.poll_interval
            check = subprocess.run(
                ["squeue", "-j", job_id, "-h"],
                capture_output=True, text=True,
            )
            if check.stdout.strip() == "":
                break  # job no longer in queue — done or failed
        else:
            return JobResult(
                status="failed",
                job_id=job_id,
                error_message=f"Job {job_id} exceeded max_wait ({slurm_cfg.max_wait}s)",
            )

        elapsed = time.time() - start

        # Check exit code via sacct
        sacct = subprocess.run(
            ["sacct", "-j", job_id, "--format=ExitCode", "--noheader", "--parsable2"],
            capture_output=True, text=True,
        )
        exit_code_str = sacct.stdout.strip().splitlines()[0] if sacct.stdout.strip() else "0:0"
        exit_code = int(exit_code_str.split(":")[0])

        stdout_text = stdout_file.read_text() if stdout_file.exists() else ""
        stderr_text = stderr_file.read_text() if stderr_file.exists() else ""

        if exit_code == 0:
            return JobResult(
                status="success",
                returncode=exit_code,
                stdout=stdout_text,
                stderr=stderr_text,
                runtime_seconds=elapsed,
                job_id=job_id,
            )
        else:
            last_err = stderr_text.strip().splitlines()[-1] if stderr_text.strip() else "Non-zero exit"
            return JobResult(
                status="failed",
                returncode=exit_code,
                stdout=stdout_text,
                stderr=stderr_text,
                runtime_seconds=elapsed,
                job_id=job_id,
                error_message=last_err,
            )


def get_executor(backend: str = "local") -> BaseExecutor:
    if backend == "local":
        return LocalExecutor()
    if backend == "slurm":
        return SlurmExecutor()
    raise ValueError(f"Unknown executor backend: {backend!r}. Supported: local, slurm")
