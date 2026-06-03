"""
Executor abstraction.
LocalExecutor runs subprocesses now.
SlurmExecutor / AwsExecutor will follow the same interface later.
"""

from __future__ import annotations

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


def get_executor(backend: str = "local") -> BaseExecutor:
    if backend == "local":
        return LocalExecutor()
    raise ValueError(f"Unknown executor backend: {backend!r}. Supported: local")
