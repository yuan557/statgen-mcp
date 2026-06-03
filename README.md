# statgen-mcp

MCP server for statistical genetics. Lets Claude orchestrate your analysis pipeline - QC, GWAS, heritability estimation, fine-mapping, and more.

## Current tools

| Tool | What it does |
|------|-------------|
| `vcf_qc` | MAF, MIND, GENO, HWE filtering via PLINK2 |
| `king_cutoff` | Relatedness filtering via PLINK2 --king-cutoff |
| `run_pca` | PCA on LD-pruned variants via PLINK2 |
| `check_tools` | Verify required binaries are installed |

## Requirements

- Python >= 3.11
- [PLINK2](https://www.cog-genomics.org/plink/2.0/) on your PATH
- [bcftools](https://samtools.github.io/bcftools/) on your PATH (required for variant source filtering)

## Installation

```bash
pip install statgen-mcp
```

Or from source:
```bash
git clone https://github.com/yuan557/statgen-mcp
cd statgen-mcp
pip install -e .
```

## Configuration

```bash
mkdir -p ~/.statgen-mcp
cp config.example.toml ~/.statgen-mcp/config.toml
```

## Claude Desktop setup

Add to `~/Library/Application Support/Claude/claude_desktop_config.json`:

```json
{
  "mcpServers": {
    "statgen-mcp": {
      "command": "statgen-mcp"
    }
  }
}
```

## Roadmap

- GWAS (REGENIE, SAIGE, PLINK2)
- Fine-mapping and conditional analysis
- Heritability (LDSC, GCTA)
- SLURM executor
- AWS Batch executor

## Compliance note

You are responsible for ensuring your use complies with applicable regulations (HIPAA, GDPR, institutional data agreements).
