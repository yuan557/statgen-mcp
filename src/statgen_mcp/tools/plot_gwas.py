"""Manhattan plot for GWAS summary statistics."""

from __future__ import annotations

from pathlib import Path

from pydantic import BaseModel, Field

from statgen_mcp.tools.plot_qq import compute_lambda_gc
from statgen_mcp.utils.validation import ValidationError, require_file


CHR_COLORS = ["#4878CF", "#85C0F9"]  # alternating per chromosome


class PlotGwasInput(BaseModel):
    summary_stats_file: str = Field(
        description="PLINK2 GWAS summary stats file (.glm.linear or .glm.logistic.hybrid)"
    )
    output_path: str = Field(description="Output PNG file path")
    trait_name: str = Field(default="Trait", description="Trait name for plot title")
    p_col: str = Field(default="P", description="Column name for p-values")
    chr_col: str = Field(default="#CHROM", description="Column name for chromosome")
    pos_col: str = Field(default="POS", description="Column name for position")
    significance_line: float = Field(default=5e-8, description="Genome-wide significance threshold")
    suggestive_line: float | None = Field(default=1e-5, description="Suggestive significance threshold")
    dpi: int = Field(default=150)
    dry_run: bool = Field(default=False)


class PlotGwasResult(BaseModel):
    status: str
    output_path: str = ""
    lambda_gc: float | None = None
    n_variants: int | None = None
    n_significant: int | None = None
    error_message: str | None = None
    command_preview: str | None = None


def run_plot_gwas(params: PlotGwasInput) -> PlotGwasResult:
    try:
        require_file(params.summary_stats_file)
    except ValidationError as e:
        return PlotGwasResult(status="failed", error_message=str(e))

    if params.dry_run:
        return PlotGwasResult(
            status="dry_run",
            command_preview=f"plot_gwas: {params.summary_stats_file} → {params.output_path}",
        )

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
        import pandas as pd
    except ImportError:
        return PlotGwasResult(
            status="failed",
            error_message="Plotting dependencies missing. Run: pip install -e '.[plots]'",
        )

    try:
        df = pd.read_csv(params.summary_stats_file, sep="\t", comment="#", low_memory=False)

        for col in [params.p_col, params.chr_col, params.pos_col]:
            if col not in df.columns:
                return PlotGwasResult(
                    status="failed",
                    error_message=f"Column '{col}' not found. Available: {list(df.columns)}",
                )

        df[params.p_col] = pd.to_numeric(df[params.p_col], errors="coerce")
        df = df.dropna(subset=[params.p_col])
        df = df[df[params.p_col] > 0]

        # Chromosome ordering: 1-22 then X/Y
        def chr_sort_key(c):
            c = str(c).replace("chr", "")
            try:
                return int(c)
            except ValueError:
                return {"X": 23, "Y": 24}.get(c, 25)

        chrs = sorted(df[params.chr_col].unique(), key=chr_sort_key)
        df["_chr_idx"] = df[params.chr_col].map({c: i for i, c in enumerate(chrs)})

        # Build x-axis positions (cumulative)
        df = df.sort_values(["_chr_idx", params.pos_col])
        chr_offsets: dict = {}
        offset = 0
        for chr_val in chrs:
            chr_offsets[chr_val] = offset
            offset += df.loc[df[params.chr_col] == chr_val, params.pos_col].max() + 1_000_000
        df["_x"] = df[params.pos_col] + df[params.chr_col].map(chr_offsets)
        df["_logp"] = -np.log10(df[params.p_col])

        lambda_gc = compute_lambda_gc(df[params.p_col])
        n_significant = int((df[params.p_col] < params.significance_line).sum())

        fig, ax = plt.subplots(figsize=(16, 5))

        for i, chr_val in enumerate(chrs):
            mask = df[params.chr_col] == chr_val
            color = CHR_COLORS[i % 2]
            ax.scatter(df.loc[mask, "_x"], df.loc[mask, "_logp"],
                       c=color, s=4, alpha=0.6, rasterized=True)

        ax.axhline(-np.log10(params.significance_line), color="red",
                   linestyle="--", linewidth=0.8, label=f"p={params.significance_line:.0e}")
        if params.suggestive_line:
            ax.axhline(-np.log10(params.suggestive_line), color="orange",
                       linestyle=":", linewidth=0.8, label=f"p={params.suggestive_line:.0e}")

        # X-axis: chromosome labels at midpoints
        chr_mids = []
        for chr_val in chrs:
            mask = df[params.chr_col] == chr_val
            mid = (df.loc[mask, "_x"].min() + df.loc[mask, "_x"].max()) / 2
            chr_mids.append((mid, str(chr_val).replace("chr", "")))
        ax.set_xticks([m for m, _ in chr_mids])
        ax.set_xticklabels([l for _, l in chr_mids], fontsize=7)

        ax.set_xlabel("Chromosome")
        ax.set_ylabel("−log₁₀(p)")
        ax.set_title(f"Manhattan Plot — {params.trait_name}  (λ_GC = {lambda_gc:.3f})")
        ax.legend(fontsize=8, loc="upper right")
        plt.tight_layout()

        Path(params.output_path).parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(params.output_path, dpi=params.dpi, bbox_inches="tight")
        plt.close()

        return PlotGwasResult(
            status="success",
            output_path=params.output_path,
            lambda_gc=round(lambda_gc, 4),
            n_variants=len(df),
            n_significant=n_significant,
        )

    except Exception as e:
        return PlotGwasResult(status="failed", error_message=str(e))
