# statgen-mcp

An MCP server that lets Claude act as an adaptive orchestrator for statistical genetics pipelines. Instead of writing a Snakefile, you tell Claude what you want — it calls the tools, reads intermediate outputs, and adjusts parameters along the way.

**Example:** "QC the imputed VCF, remove related samples, run PCA, then run GWAS on BMI adjusting for age, sex, and PCs."
Claude will chain `vcf_qc → king_cutoff → filter_samples → run_pca → run_gwas`, pick `n_pcs` from the eigenvalue curve, and flag anything unexpected in the QC logs.

## Tools

| Tool | What it does |
|------|-------------|
| `vcf_qc` | MAF, MIND, GENO, HWE filtering via PLINK2. Pre-filters VCF by variant source (genotyped / imputed INFO score / WGS) via bcftools. |
| `king_cutoff` | Relatedness analysis via PLINK2 KING. Report-only mode or filter mode (removes one sample per related pair). |
| `het_check` | Heterozygosity QC via PLINK2 `--het`. Flags contaminated samples (excess het) and inbred samples (deficient het) by F statistic SD cutoff. Writes a remove file for `filter_samples`. |
| `filter_samples` | Keep or remove samples by FID/IID list — ancestry subsets, PCA outliers, kinship exclusions. |
| `run_pca` | LD-prune then PCA via PLINK2. Optionally excludes high-LD regions (MHC, inversions). Outputs `.eigenvec` / `.eigenval`. |
| `run_gwas` | PLINK2 `--glm` per chromosome. Merges user covariates + PCA eigenvec internally. Firth fallback for binary traits. |
| `plot_pca` | PCA scatter plots coloured by population label. |
| `plot_gwas` | Manhattan plot with genome-wide / suggestive significance lines and lambda GC. |
| `plot_qq` | QQ plot with lambda GC and lambda_1000 annotation. |
| `check_tools` | Verify required binaries (plink2, bcftools) are installed and show versions. |

## Intended pipeline order

```
vcf_qc (genotyped or imputed INFO >= 0.99)
  het_check (remove contaminated / inbred samples)
  king_cutoff (remove cryptic relatives)
  filter_samples (ancestry subset, outlier removal)
  run_pca
  run_gwas (imputed INFO >= 0.3–0.8, per chromosome)
```

## Requirements

- Python >= 3.11
- [PLINK2](https://www.cog-genomics.org/plink/2.0/) on your PATH
- [bcftools](https://samtools.github.io/bcftools/) on your PATH (required for variant source filtering)

## Installation

```bash
pip install statgen-mcp
```

With plotting support:
```bash
pip install "statgen-mcp[plots]"
```

Or from source:
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

Output files go to `~/statgen-mcp-runs/` by default (one timestamped subdirectory per tool call). SLURM submission is supported — see `config.example.toml`.

## Claude Desktop / Claude Code setup

Add to your MCP config:

```json
{
  "mcpServers": {
    "statgen-mcp": {
      "command": "statgen-mcp"
    }
  }
}
```

For Claude Code (CLI), add to `.mcp.json` in your project root or `~/.claude.json` globally.

## Roadmap

- `merge_plink` — merge per-chromosome PLINK files for genome-wide PCA/KING
- `sanity_check` — sex concordance across phenotype file, .fam, and chrX
- `plot_kinship` — kinship histogram and IBS0 vs kinship scatter from `.kin0`
- PRS (clump + threshold, internal and external summary stats)
- Fine-mapping and conditional analysis
- Heritability (LDSC, GCTA)
- AWS Batch executor

## Compliance note

You are responsible for ensuring your use complies with applicable regulations (HIPAA, GDPR, institutional data agreements).
