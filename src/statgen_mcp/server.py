"""
statgen-mcp: MCP server for statistical genetics.
Entry point - registers tools with the MCP SDK.
"""

from __future__ import annotations

import json

from mcp.server import Server
from mcp.server.stdio import stdio_server
from mcp.types import TextContent, Tool

from statgen_mcp.tools.filter_samples import FilterSamplesInput, run_filter_samples
from statgen_mcp.tools.gwas import GwasInput, run_gwas
from statgen_mcp.tools.het_check import HetCheckInput, run_het_check
from statgen_mcp.tools.plot_gwas import PlotGwasInput, run_plot_gwas
from statgen_mcp.tools.plot_pca import PlotPcaInput, run_plot_pca
from statgen_mcp.tools.plot_qq import PlotQqInput, run_plot_qq
from statgen_mcp.tools.king import KingInput, run_king_cutoff
from statgen_mcp.tools.pca import PcaInput, run_pca
from statgen_mcp.tools.data_qc import DataQcInput, run_data_qc
from statgen_mcp.tools.merge_plink import MergePlinkInput, run_merge_plink

app = Server("statgen-mcp")


@app.list_tools()
async def list_tools() -> list[Tool]:
    return [
        Tool(
            name="merge_plink",
            description=(
                "Merge multiple PLINK datasets into one using PLINK2 --pmerge-list. "
                "Typical use: combine per-chromosome PLINK files into a single genome-wide dataset "
                "before running data_qc, king_cutoff, or run_pca. "
                "Accepts PLINK1 (bed/bim/fam) or PLINK2 native (pgen/psam/pvar) input; "
                "output format matches input format. At least two prefixes required."
            ),
            inputSchema=MergePlinkInput.model_json_schema(),
        ),
        Tool(
            name="data_qc",
            description=(
                "Filter a VCF or PLINK dataset for quality control. "
                "Accepts VCF/VCF.gz, PLINK1 (bed/bim/fam), or PLINK2 native (pgen/psam/pvar) input. "
                "Applies MAF, per-sample missingness (MIND), per-variant missingness (GENO), "
                "and Hardy-Weinberg equilibrium (HWE) filters using PLINK2. "
                "For VCF input, optionally pre-filters variants by source: "
                "'genotyped' (directly typed array variants only), "
                "'imputed' (imputation INFO score >= threshold), or 'wgs' (no INFO filter). "
                "Outputs a clean PLINK bed/bim/fam dataset."
            ),
            inputSchema=DataQcInput.model_json_schema(),
        ),
        Tool(
            name="king_cutoff",
            description=(
                "Relatedness analysis using PLINK2 KING. "
                "By default (report_only=False), removes one sample from each pair above the kinship cutoff. "
                "With report_only=True, reports all related pairs above the cutoff without removing anyone. "
                "Typical GWAS cutoff: 0.0442 (3rd-degree relatives). "
                "Input must be a PLINK prefix. "
                "For imputed input, prepare with data_qc using variant_source='imputed' and min_info_score >= 0.99, "
                "or use variant_source='genotyped' to restrict to directly typed variants."
            ),
            inputSchema=KingInput.model_json_schema(),
        ),
        Tool(
            name="run_pca",
            description=(
                "Run principal component analysis on a PLINK dataset using PLINK2. "
                "LD-prunes first by default. "
                "Outputs .eigenvec and .eigenval files. "
                "Use PCs as covariates in downstream GWAS. "
                "Input should be genotyped variants or imputed variants with INFO >= 0.99 — "
                "prepare with data_qc using variant_source='genotyped' or variant_source='imputed' min_info_score=0.99."
            ),
            inputSchema=PcaInput.model_json_schema(),
        ),
        Tool(
            name="plot_pca",
            description=(
                "Plot PCA results as scatter plots (PC1 vs PC2, PC3 vs PC4 by default). "
                "Optionally colour points by population label from a separate file. "
                "Input is the .eigenvec file from run_pca."
            ),
            inputSchema=PlotPcaInput.model_json_schema(),
        ),
        Tool(
            name="plot_gwas",
            description=(
                "Manhattan plot for GWAS summary statistics. "
                "Draws genome-wide significance (5e-8) and suggestive (1e-5) lines. "
                "Annotates plot with lambda GC. Returns n_significant variants. "
                "Input is a PLINK2 .glm.linear or .glm.logistic.hybrid summary stats file."
            ),
            inputSchema=PlotGwasInput.model_json_schema(),
        ),
        Tool(
            name="plot_qq",
            description=(
                "QQ plot for GWAS summary statistics. "
                "Computes and annotates lambda GC (and lambda_1000 for binary traits). "
                "Input is a PLINK2 .glm.linear or .glm.logistic.hybrid summary stats file."
            ),
            inputSchema=PlotQqInput.model_json_schema(),
        ),
        Tool(
            name="run_gwas",
            description=(
                "Run genome-wide association study using PLINK2 --glm. "
                "Runs on the full dataset in a single call — PLINK2 outputs one summary stats file per phenotype covering all chromosomes. "
                "Merges user covariates and PCA eigenvec into a combined covariate file. "
                "Auto-detects binary vs quantitative trait; use firth_fallback=True (default) for binary traits. "
                "Input should be prepared with data_qc. PCs from run_pca should use INFO >= 0.99 or genotyped variants. "
                "For GWAS association input, use a separate data_qc run with a lower INFO threshold (e.g. 0.3–0.8)."
            ),
            inputSchema=GwasInput.model_json_schema(),
        ),
        Tool(
            name="filter_samples",
            description=(
                "Filter samples from a PLINK dataset using --keep and/or --remove. "
                "Accepts any FID/IID sample list: ancestry subsets (e.g. EUR.txt), "
                "PCA outliers, kinship exclusions (pass .king.cutoff.out.id directly), "
                "or user-specified lists. At least one of keep_file or remove_file is required."
            ),
            inputSchema=FilterSamplesInput.model_json_schema(),
        ),
        Tool(
            name="het_check",
            description=(
                "Heterozygosity QC via PLINK2 --het. "
                "Flags samples deviating more than sd_threshold SD from the mean F statistic. "
                "Excess heterozygosity (F << mean) indicates sample contamination. "
                "Deficient heterozygosity (F >> mean) indicates inbreeding or population structure. "
                "Run on genotyped variants after data_qc, before king_cutoff. "
                "Optionally writes a remove file ready for filter_samples."
            ),
            inputSchema=HetCheckInput.model_json_schema(),
        ),
        Tool(
            name="check_tools",
            description=(
                "Check which required external tools (plink2, bcftools) are installed "
                "and available on PATH. Returns version info where available."
            ),
            inputSchema={"type": "object", "properties": {}, "required": []},
        ),
    ]


@app.call_tool()
async def call_tool(name: str, arguments: dict) -> list[TextContent]:
    try:
        if name == "merge_plink":
            params = MergePlinkInput(**arguments)
            result = run_merge_plink(params)
            return [TextContent(type="text", text=result.model_dump_json(indent=2))]

        elif name == "data_qc":
            params = DataQcInput(**arguments)
            result = run_data_qc(params)
            return [TextContent(type="text", text=result.model_dump_json(indent=2))]

        elif name == "plot_pca":
            params = PlotPcaInput(**arguments)
            result = run_plot_pca(params)
            return [TextContent(type="text", text=result.model_dump_json(indent=2))]

        elif name == "plot_gwas":
            params = PlotGwasInput(**arguments)
            result = run_plot_gwas(params)
            return [TextContent(type="text", text=result.model_dump_json(indent=2))]

        elif name == "plot_qq":
            params = PlotQqInput(**arguments)
            result = run_plot_qq(params)
            return [TextContent(type="text", text=result.model_dump_json(indent=2))]

        elif name == "run_gwas":
            params = GwasInput(**arguments)
            result = run_gwas(params)
            return [TextContent(type="text", text=result.model_dump_json(indent=2))]

        elif name == "filter_samples":
            params = FilterSamplesInput(**arguments)
            result = run_filter_samples(params)
            return [TextContent(type="text", text=result.model_dump_json(indent=2))]

        elif name == "het_check":
            params = HetCheckInput(**arguments)
            result = run_het_check(params)
            return [TextContent(type="text", text=result.model_dump_json(indent=2))]

        elif name == "king_cutoff":
            params = KingInput(**arguments)
            result = run_king_cutoff(params)
            return [TextContent(type="text", text=result.model_dump_json(indent=2))]

        elif name == "run_pca":
            params = PcaInput(**arguments)
            result = run_pca(params)
            return [TextContent(type="text", text=result.model_dump_json(indent=2))]

        elif name == "check_tools":
            return [TextContent(type="text", text=json.dumps(_check_tools(), indent=2))]

        else:
            return [TextContent(type="text", text=json.dumps({"error": f"Unknown tool: {name}"}))]

    except Exception as e:
        return [TextContent(type="text", text=json.dumps({"error": str(e)}))]


def _check_tools() -> dict:
    import shutil
    import subprocess

    tools = ["plink2", "bcftools"]
    results = {}

    for tool in tools:
        path = shutil.which(tool)
        if path is None:
            results[tool] = {"available": False, "path": None, "version": None}
        else:
            version = None
            try:
                out = subprocess.run([tool, "--version"], capture_output=True, text=True, timeout=5)
                version = (out.stdout or out.stderr).strip().splitlines()[0]
            except Exception:
                pass
            results[tool] = {"available": True, "path": path, "version": version}

    return results


async def _main() -> None:
    async with stdio_server() as (read_stream, write_stream):
        await app.run(read_stream, write_stream, app.create_initialization_options())


def main() -> None:
    import asyncio
    asyncio.run(_main())


if __name__ == "__main__":
    main()
