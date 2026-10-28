# statgen-mcp

[![Tests](https://github.com/yuan557/statgen-mcp/actions/workflows/tests.yml/badge.svg)](https://github.com/yuan557/statgen-mcp/actions/workflows/tests.yml)
[![Integration Tests](https://github.com/yuan557/statgen-mcp/actions/workflows/integration.yml/badge.svg)](https://github.com/yuan557/statgen-mcp/actions/workflows/integration.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

An MCP server that exposes PLINK2 and bcftools as tools so an AI agent can orchestrate GWAS pipelines in plain English - no Snakefile required.

```
> QC the imputed VCF, remove related samples, run PCA, then run GWAS on BMI adjusting for age, sex, and PCs.
```

The agent chains `data_qc → king_cutoff → filter_samples → run_pca → run_gwas`, reads intermediate outputs, and adjusts parameters along the way.

![Pipeline output panel](docs/images/pipeline_panel.png)

*504 East Asian samples, 1.23M variants, simulated binary trait. λ_GC = 1.038, 4/5 causal SNPs recovered.*

## Tools

| Tool | What it does |
|------|-------------|
| `merge_plink` | Merge per-chromosome PLINK datasets |
| `data_qc` | MAF / MIND / GENO / HWE filtering. Accepts VCF, PLINK1, PLINK2 |
| `king_cutoff` | Relatedness filtering or report-only mode |
| `het_check` | Flag contamination outliers by heterozygosity F statistic |
| `filter_samples` | Keep or remove samples by FID/IID list |
| `run_pca` | LD-prune then PCA |
| `run_gwas` | PLINK2 `--glm` with Firth fallback for binary traits |
| `plot_pca` | PCA scatter coloured by population label |
| `plot_gwas` | Manhattan plot with λ_GC |
| `plot_qq` | QQ plot with λ_GC annotation |

All tools accept PLINK1 (`.bed/.bim/.fam`) and PLINK2 native (`.pgen/.psam/.pvar`). Format is auto-detected; output matches input.

## Installation

**Conda (local / workstation):**
```bash
conda create -n statgen-mcp -c bioconda -c conda-forge plink2 bcftools python>=3.11
conda activate statgen-mcp
pip install -e ".[plots]"
```

**HPC (module + venv):**
```bash
module load plink2 bcftools
python -m venv ~/.venv/statgen-mcp && source ~/.venv/statgen-mcp/bin/activate
pip install -e ".[plots]"
```

## Agent setup

### Claude Code

Add to `.mcp.json` in your project root:
```json
{
  "mcpServers": {
    "statgen-mcp": { "command": "statgen-mcp" }
  }
}
```
Run `claude` in your project directory and describe the analysis in plain English.

### AWS EC2 (self-hosted Qwen via vLLM)

Run Qwen3-4B-AWQ on a `g4dn.xlarge` spot instance (~$0.16/hr). Data stays local - only LLM calls go to EC2.

On EC2 (Deep Learning AMI):
```bash
python -m venv vllm_env && source vllm_env/bin/activate
pip install -r deploy/requirements.txt
bash deploy/start_vllm.sh
```

Add to `~/.statgen-mcp/config.toml`:
```toml
[llm]
base_url = "http://<ec2-ip>:8000/v1"
model    = "Qwen/Qwen3-4B-AWQ"
api_key  = "none"
```

Then run:
```bash
pip install -e ".[llm]"
python deploy/run_pipeline.py "QC /data/cohort variant_source=genotyped, then run PCA"
```

### AWS EC2 (CPU - no GPU required)

Run Qwen3-4B on a `c7i.2xlarge` (8 vCPU, 16GB RAM, ~$0.34/hr) using llama.cpp. Expect ~5-10 tok/s - slow for chat, fine for driving a pipeline where tool calls take minutes.

On EC2 (Amazon Linux 2023 AMI):
```bash
python -m venv llama_env && source llama_env/bin/activate
pip install -r deploy/requirements-cpu.txt
bash deploy/start_llama_cpp.sh   # downloads Qwen3-1.7B GGUF on first run (~1GB)
```

Add to `~/.statgen-mcp/config.toml`:
```toml
[llm]
base_url = "http://<ec2-ip>:8000/v1"
model    = "Qwen/Qwen3-4B-GGUF"
api_key  = "none"
```

Then run `deploy/run_pipeline.py` as above - the API is identical to the GPU setup.

### Local GPU

Same as the GPU EC2 option - run `deploy/start_vllm.sh` on any Linux machine with an NVIDIA GPU (≥ 6GB VRAM), point `base_url` at `http://localhost:8000/v1`.

## Configuration

```bash
mkdir -p ~/.statgen-mcp && cp config.example.toml ~/.statgen-mcp/config.toml
```

Key settings: `output_dir`, `[executor] backend` (`local` or `slurm`), binary paths under `[tools]`. See `config.example.toml` for SLURM and LLM options.

## Acknowledgements

Pipeline validation uses the 1000 Genomes EAS dataset from [GWASTutorial](https://github.com/Cloufield/GWASTutorial) by Cloufield et al.

## Compliance note

You are responsible for ensuring your use complies with applicable regulations (HIPAA, GDPR, institutional data agreements).
