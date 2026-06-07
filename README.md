# statgen-mcp

An MCP server that wraps **PLINK2** and **bcftools** and lets any agentic AI (e.g. Claude, Gemini, Codex) act as an adaptive orchestrator for statistical genetics pipelines. Instead of writing a Snakefile, you tell the agent what you want — it calls the tools, reads intermediate outputs, and adjusts parameters along the way.

**Example:** "QC the imputed VCF, remove related samples, run PCA, then run GWAS on BMI adjusting for age, sex, and PCs."
The agent will chain `vcf_qc`, `king_cutoff`, `filter_samples`, `run_pca`, and `run_gwas`, pick `n_pcs` from the eigenvalue curve, and flag anything unexpected in the QC logs.

## Implemented tools

| Tool | What it does |
|------|-------------|
| `vcf_qc` | MAF, MIND, GENO, HWE filtering via PLINK2. Pre-filters VCF by variant source (genotyped / imputed INFO score / WGS) via bcftools. |
| `king_cutoff` | Relatedness analysis via PLINK2 KING. Report-only mode or filter mode (removes one sample per related pair). |
| `het_check` | Heterozygosity QC via PLINK2 `--het`. Flags contaminated samples (excess het) and inbred samples (deficient het) by F statistic SD cutoff. Writes a remove file for `filter_samples`. |
| `filter_samples` | Keep or remove samples by FID/IID list: ancestry subsets, PCA outliers, kinship exclusions. |
| `run_pca` | LD-prune then PCA via PLINK2. Optionally excludes high-LD regions (MHC, inversions). Outputs `.eigenvec` / `.eigenval`. |
| `run_gwas` | PLINK2 `--glm` on the full dataset. Merges user covariates + PCA eigenvec internally. Firth fallback for binary traits. |
| `plot_pca` | PCA scatter plots coloured by population label. |
| `plot_gwas` | Manhattan plot with genome-wide / suggestive significance lines and lambda GC. |
| `plot_qq` | QQ plot with lambda GC annotation. |
| `check_tools` | Verify required binaries (plink2, bcftools) are installed and show versions. |

All tools accept both **PLINK1** (`.bed/.bim/.fam`) and **PLINK2 native** (`.pgen/.psam/.pvar`) format. Format is auto-detected from the input prefix; output format matches input format. `vcf_qc` additionally accepts VCF/VCF.gz and always outputs PLINK1 bed format.

## Intended pipeline order

**Step 1: Variant QC**
Run `vcf_qc` on genotyped or high-confidence imputed variants (INFO >= 0.99). Outputs a clean PLINK dataset.

**Step 2: Sample QC** (run all, collect remove files)
- `het_check`: flag contaminated samples (excess heterozygosity) and inbred samples (deficient heterozygosity)
- `king_cutoff`: identify cryptic relatives
- *(optional)* manual ancestry filtering using population labels

**Step 3: Apply exclusions**
Pass the remove files from step 2 to `filter_samples`.

**Step 4: PCA**
Run `run_pca` on the clean, filtered dataset. Review `plot_pca`. Choose n_pcs from the eigenvalue curve.

**Step 5: GWAS**
Run `run_gwas` on imputed variants (INFO >= 0.3-0.8). Pass PCs as covariates. Review `plot_gwas` and `plot_qq`; check lambda GC.

## Requirements

- Python >= 3.11
- [PLINK2](https://www.cog-genomics.org/plink/2.0/) on your PATH
- [bcftools](https://samtools.github.io/bcftools/) on your PATH (required for variant source filtering)

## Installation

### Local / workstation (conda recommended)

conda can install the Python package and all binary dependencies together via [bioconda](https://bioconda.github.io/):

```bash
conda create -n statgen-mcp -c bioconda -c conda-forge plink2 bcftools python>=3.11
conda activate statgen-mcp
pip install "statgen-mcp[plots]"
```

### HPC clusters (e.g. Compute Canada / Alliance)

Many HPC clusters discourage conda due to Lustre filesystem overhead. Use the cluster module system for binaries and a Python virtualenv for the package:

```bash
module load plink2 bcftools   # exact module names vary by cluster
python -m venv ~/.venv/statgen-mcp
source ~/.venv/statgen-mcp/bin/activate
pip install "statgen-mcp[plots]"
```

Add the `module load` and `source` lines to your `~/.bashrc` (or job script) so the environment is active when the MCP server starts.

### From source

```bash
git clone https://github.com/yuan557/statgen-mcp
cd statgen-mcp
pip install -e ".[plots]"
```

## Configuration

```bash
mkdir -p ~/.statgen-mcp
cp config.example.toml ~/.statgen-mcp/config.toml
```

Output files go to `~/statgen-mcp-runs/` by default (one timestamped subdirectory per tool call). SLURM submission is supported - see `config.example.toml`.

## Agent setup

Works with any MCP-compatible agent. Examples:

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

**Claude Desktop** - add the same block to `~/Library/Application Support/Claude/claude_desktop_config.json`.

Other agents (Gemini CLI, Codex CLI, etc.) - consult their MCP client documentation for the equivalent config.

## Roadmap (not yet implemented)

- `merge_plink`: merge per-chromosome PLINK files for genome-wide PCA/KING
- `sanity_check`: sex concordance across phenotype file, .fam, and chrX
- `plot_kinship`: kinship histogram and IBS0 vs kinship scatter from `.kin0`
- PRS (clump + threshold, internal and external summary stats)
- Fine-mapping and conditional analysis
- Heritability (LDSC, GCTA)
- AWS Batch executor

## Compliance note

You are responsible for ensuring your use complies with applicable regulations (HIPAA, GDPR, institutional data agreements).
