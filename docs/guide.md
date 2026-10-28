# statgen-mcp guide

Full technical reference. For the short pitch see the [README](../README.md).

## Contents

- [Tools](#tools)
- [Input formats](#input-formats)
- [Installation](#installation)
- [Configuration](#configuration)
- [Agent setup](#agent-setup)
- [Pipeline order](#pipeline-order)
- [Validation](#validation)

## Tools

| Tool | What it does |
|------|-------------|
| `merge_plink` | Merge per-chromosome PLINK datasets into one genome-wide file via PLINK2 `--pmerge-list`. Accepts PLINK1 or PLINK2 native input. |
| `data_qc` | MAF, MIND, GENO, HWE filtering via PLINK2. Accepts VCF, PLINK1, or PLINK2 native input. Pre-filters VCF by variant source (genotyped / imputed INFO score / WGS) via bcftools. |
| `king_cutoff` | Relatedness analysis via PLINK2 KING. Report-only mode or filter mode (removes one sample per related pair, writes `.king.cutoff.out.id` for `filter_samples`). |
| `het_check` | Heterozygosity QC via PLINK2 `--het`. Flags outliers that indicate contamination by F statistic SD cutoff (default ±3 SD). Writes a remove file for `filter_samples`. |
| `filter_samples` | Keep or remove samples by FID/IID list: ancestry subsets, PCA outliers, kinship exclusions, het outliers. |
| `run_pca` | LD-prune then PCA via PLINK2. Optionally excludes high-LD regions (MHC, inversions). Outputs `.eigenvec` / `.eigenval`. |
| `run_gwas` | PLINK2 `--glm` on the full dataset. Merges user covariates + PCA eigenvec internally. Firth fallback for binary traits. |
| `plot_pca` | PCA scatter plots coloured by population label. |
| `plot_gwas` | Manhattan plot with genome-wide / suggestive significance lines and lambda GC. |
| `plot_qq` | QQ plot with lambda GC annotation. |
| `check_tools` | Verify required binaries (plink2, bcftools) are installed and show versions. |

## Input formats

All tools accept both **PLINK1** (`.bed/.bim/.fam`) and **PLINK2 native** (`.pgen/.psam/.pvar`) format. Format is auto-detected from the input prefix; output format matches input format. `data_qc` additionally accepts VCF/VCF.gz input, which always outputs PLINK1 bed format.

## Installation

### Local / workstation (conda recommended)

conda can install the Python package and all binary dependencies together via [bioconda](https://bioconda.github.io/):

```bash
conda create -n statgen-mcp -c bioconda -c conda-forge plink2 bcftools python>=3.11
conda activate statgen-mcp
git clone https://github.com/yuan557/statgen-mcp
cd statgen-mcp
pip install -e ".[plots]"
```

### HPC clusters (e.g. Compute Canada / Alliance)

Many HPC clusters discourage conda due to Lustre filesystem overhead. Use the cluster module system for binaries and a Python virtualenv for the package:

```bash
module load plink2 bcftools   # exact module names vary by cluster
git clone https://github.com/yuan557/statgen-mcp
cd statgen-mcp
python -m venv ~/.venv/statgen-mcp
source ~/.venv/statgen-mcp/bin/activate
pip install -e ".[plots]"
```

Add the `module load` and `source` lines to your `~/.bashrc` (or job script) so the environment is active when the MCP server starts.

## Configuration

```bash
mkdir -p ~/.statgen-mcp
cp config.example.toml ~/.statgen-mcp/config.toml
```

The key settings to review:

```toml
# Where outputs are written (one timestamped subdirectory per tool call)
output_dir = "~/statgen-mcp-runs"

[executor]
backend = "local"   # or "slurm" on HPC

[tools]
plink2 = "plink2"       # full path if not on PATH, e.g. "/usr/local/bin/plink2"
bcftools = "bcftools"   # required for VCF variant source filtering
```

SLURM options (partition, account, memory, time) are in `config.example.toml`.

## Agent setup

Works with any MCP-compatible agent.

**Claude Code (CLI)** - add to `.mcp.json` in your project root:
```json
{
  "mcpServers": {
    "statgen-mcp": {
      "command": "statgen-mcp"
    }
  }
}
```

**Claude Desktop** - add the same block to the config file for your OS:
- macOS: `~/Library/Application Support/Claude/claude_desktop_config.json`
- Windows: `%APPDATA%\Claude\claude_desktop_config.json`
- Linux: `~/.config/Claude/claude_desktop_config.json`

Other agents (Gemini CLI, Codex CLI, etc.) - consult their MCP client documentation for the equivalent config.

### Self-hosted LLM backend (experimental)

Instead of Claude, you can drive the pipeline with a self-hosted open-source model exposing an OpenAI-compatible API. Configure it in `~/.statgen-mcp/config.toml`:

```toml
[llm]
base_url = "http://<host>:8000/v1"
model    = "Qwen/Qwen3-4B-AWQ"
api_key  = "none"
```

Then run `deploy/run_pipeline.py` (install with `pip install -e ".[llm]"`). Data stays local; only LLM calls go to the host.

Status:

- **GPU (vLLM)** - in progress. A `g4dn.xlarge` (T4, 16 GB VRAM) running Qwen3-4B-AWQ is the intended path.
- **CPU (llama.cpp)** - tested and shelved. An 8 GB instance is not enough (llama.cpp pre-allocates ~7.5 GB, leaving no headroom and getting OOM-killed), and Qwen3-1.7B is too small to orchestrate a multi-step pipeline reliably - it loops on errors instead of adjusting. Needs >8 GB RAM and a larger model (Qwen3-4B+) to be usable.

## Usage example

With `.mcp.json` in place, start Claude Code in your project directory and describe what you want in plain English:

```
$ claude

> I have an imputed VCF at /data/cohort.vcf.gz and a phenotype file at /data/pheno.txt
  (columns: FID IID BMI AGE SEX). QC the imputed variants (INFO >= 0.3), remove
  related samples and het outliers, run PCA, then run GWAS on BMI adjusting for
  age, sex, and PCs.
```

Claude calls the tools in sequence - `data_qc`, `het_check`, `king_cutoff`, `filter_samples`, `run_pca`, `run_gwas`, `plot_gwas`, `plot_qq` - reading intermediate outputs at each step. It picks the number of PCs from the eigenvalue curve, checks λ_GC after GWAS, and flags anything unexpected in the QC logs. You supervise; it executes.

Example summary at the end:

```
Pipeline complete.
- Variants: 8.2M → 6.1M after QC (INFO >= 0.3, MAF 1%, geno 5%, HWE 1e-6)
- Samples: 4,892 → 4,801 (45 het outliers, 46 related removed)
- PCA: 10 PCs explain 87% of variance
- GWAS: λ_GC = 1.02, 3 genome-wide significant loci
  Plots: /data/results/.../manhattan.png, qq.png
```

## Pipeline order

### Step 1: Variant QC

Run `data_qc` on genotyped or high-confidence imputed variants (INFO >= 0.99). Outputs a clean PLINK dataset used for PCA and kinship.

### Step 2: Sample QC (run all, collect remove files)

- `het_check` - flag outliers that indicate contamination
- `king_cutoff` - identify cryptic relatives
- *(optional)* manual ancestry filtering using population labels

### Step 3: Apply exclusions

Pass the remove files from step 2 to `filter_samples`:

- `.king.cutoff.out.id` from `king_cutoff` to remove related individuals
- `het_outliers_remove.txt` from `het_check` to remove het outliers
- *(optional)* ancestry subset file

### Step 4: PCA

Run `run_pca` on the clean, filtered dataset. Review `plot_pca`. Choose `n_pcs` from the eigenvalue curve.

### Step 5: GWAS

The Step 1 dataset is for PCA and kinship only. Run a second `data_qc` to produce the GWAS input - the right settings depend on your data type:

- **Imputed:** use a lower INFO threshold (e.g. `min_info_score=0.3`, `variant_source="imputed"`). The 0.99 dataset discards most imputed variants.
- **WGS:** `variant_source="wgs"` - no INFO filtering. Adjust MAF and other filters as appropriate.
- **Genotyped array:** the Step 1 dataset is usually fine as-is; no second run needed.

Then run `run_gwas` on this dataset, passing PCs as covariates. Review `plot_gwas` and `plot_qq`; check λ_GC.

## Validation

The README example is a run on 504 East Asian samples from the 1000 Genomes Project using a simulated binary trait with 5 planted causal SNPs (effect size OR=3, MAF > 1%).

Planted causal SNPs are real genotyped variants with a synthetic phenotype constructed so that exactly those 5 SNPs drive case-control status. Because the ground truth is known, you can check whether the GWAS recovers them - a stronger validation than λ_GC alone, since it tests whether the pipeline produces correct hits, not just well-calibrated p-values.

Results:

- Variants: 1.23M → 1.03M after QC (MAF 1%, missingness, HWE)
- Samples: 504 → 495 after kinship filtering (7 related pairs)
- PCA: clean EAS subpopulation separation (CDX/CHB/CHS/JPT/KHV)
- GWAS: λ_GC = 1.038, 57 genome-wide hits, 4/5 causal SNPs recovered
- The missed SNP (chr8) had MAF = 0.2% and was correctly removed by the MAF filter
