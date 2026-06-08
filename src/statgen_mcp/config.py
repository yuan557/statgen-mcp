"""
Configuration for statgen-mcp.
Reads from ~/.statgen-mcp/config.toml or STATGEN_MCP_CONFIG env var.
"""

from __future__ import annotations

import os
import tomllib
from dataclasses import dataclass, field
from pathlib import Path


DEFAULT_CONFIG_PATH = Path.home() / ".statgen-mcp" / "config.toml"


@dataclass
class ExecutorConfig:
    backend: str = "local"  # local | slurm | aws | gcp | azure


@dataclass
class SlurmConfig:
    partition: str = "normal"
    time: str = "04:00:00"
    default_memory: str = "16G"
    default_cpus: int = 4
    account: str = ""
    poll_interval: int = 30    # seconds between squeue checks
    max_wait: int = 86400      # seconds before giving up (24h)


@dataclass
class AwsConfig:
    region: str = "us-east-1"
    job_queue: str = ""
    job_definition: str = ""


@dataclass
class ToolsConfig:
    plink2: str = "plink2"
    bcftools: str = "bcftools"


@dataclass
class Config:
    executor: ExecutorConfig = field(default_factory=ExecutorConfig)
    slurm: SlurmConfig = field(default_factory=SlurmConfig)
    aws: AwsConfig = field(default_factory=AwsConfig)
    tools: ToolsConfig = field(default_factory=ToolsConfig)
    output_dir: str = str(Path.home() / "statgen-mcp-runs")


def load_config() -> Config:
    config_path = Path(os.environ.get("STATGEN_MCP_CONFIG", DEFAULT_CONFIG_PATH))

    if not config_path.exists():
        return Config()

    with open(config_path, "rb") as f:
        raw = tomllib.load(f)

    cfg = Config()

    if "executor" in raw:
        cfg.executor = ExecutorConfig(**raw["executor"])
    if "slurm" in raw:
        cfg.slurm = SlurmConfig(**raw["slurm"])
    if "aws" in raw:
        cfg.aws = AwsConfig(**raw["aws"])
    if "tools" in raw:
        cfg.tools = ToolsConfig(**raw["tools"])
    if "output_dir" in raw:
        cfg.output_dir = raw["output_dir"]

    return cfg


_config: Config | None = None


def get_config() -> Config:
    global _config
    if _config is None:
        _config = load_config()
    return _config
