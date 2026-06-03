"""
statgen-mcp: MCP server for statistical genetics.
Entry point - registers tools with the MCP SDK.
"""

from __future__ import annotations

import json

from mcp.server import Server
from mcp.server.stdio import stdio_server
from mcp.types import TextContent, Tool

from statgen_mcp.tools.king import KingInput, run_king_cutoff
from statgen_mcp.tools.pca import PcaInput, run_pca
from statgen_mcp.tools.vcf_qc import VcfQcInput, run_vcf_qc

app = Server("statgen-mcp")


@app.list_tools()
async def list_tools() -> list[Tool]:
    return [
        Tool(
            name="vcf_qc",
            description=(
                "Filter a VCF or PLINK dataset for quality control. "
                "Applies MAF, per-sample missingness (MIND), per-variant missingness (GENO), "
                "and Hardy-Weinberg equilibrium (HWE) filters using PLINK2. "
                "For VCF input, optionally pre-filters variants by source: "
                "'genotyped' (directly typed array variants only), "
                "'imputed' (imputation INFO score >= threshold), or 'wgs' (no INFO filter). "
                "Outputs a clean PLINK bed/bim/fam dataset."
            ),
            inputSchema=VcfQcInput.model_json_schema(),
        ),
        Tool(
            name="king_cutoff",
            description=(
                "Remove related individuals using PLINK2 --king-cutoff. "
                "One sample from each related pair above the kinship cutoff is removed. "
                "Typical GWAS cutoff: 0.0442 (3rd-degree relatives). "
                "Input must be a PLINK prefix."
            ),
            inputSchema=KingInput.model_json_schema(),
        ),
        Tool(
            name="run_pca",
            description=(
                "Run principal component analysis on a PLINK dataset using PLINK2. "
                "LD-prunes first by default. "
                "Outputs .eigenvec and .eigenval files. "
                "Use PCs as covariates in downstream GWAS."
            ),
            inputSchema=PcaInput.model_json_schema(),
        ),
        Tool(
            name="check_tools",
            description=(
                "Check which required external tools (plink2, flashpca) are installed "
                "and available on PATH. Returns version info where available."
            ),
            inputSchema={"type": "object", "properties": {}, "required": []},
        ),
    ]


@app.call_tool()
async def call_tool(name: str, arguments: dict) -> list[TextContent]:
    try:
        if name == "vcf_qc":
            params = VcfQcInput(**arguments)
            result = run_vcf_qc(params)
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

    tools = ["plink2", "bcftools", "flashpca"]
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
