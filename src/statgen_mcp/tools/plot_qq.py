"""QQ plot with lambda GC for GWAS summary statistics."""

from __future__ import annotations

from pathlib import Path

from pydantic import BaseModel, Field

from statgen_mcp.utils.validation import ValidationError, require_file


class PlotQqInput(BaseModel):
    summary_stats_file: str = Field(
        description="PLINK2 GWAS summary stats file (.glm.linear or .glm.logistic.hybrid)"
    )
    output_path: str = Field(description="Output PNG file path")
    p_col: str = Field(default="P", description="Column name for p-values")
    dpi: int = Field(default=150)
    dry_run: bool = Field(default=False)


class PlotQqResult(BaseModel):
    status: str
    output_path: str = ""
    lambda_gc: float | None = None
    n_variants: int | None = None
    error_message: str | None = None
    command_preview: str | None = None


def compute_lambda_gc(pvalues) -> float:
    import numpy as np
    from scipy.stats import chi2
    pvalues = pvalues.dropna()
    pvalues = pvalues[pvalues > 0]
    chisq = chi2.ppf(1 - pvalues, df=1)
    return float(np.median(chisq) / chi2.ppf(0.5, df=1))


def run_plot_qq(params: PlotQqInput) -> PlotQqResult:
    try:
        require_file(params.summary_stats_file)
    except ValidationError as e:
        return PlotQqResult(status="failed", error_message=str(e))

    if params.dry_run:
        return PlotQqResult(
            status="dry_run",
            command_preview=f"plot_qq: {params.summary_stats_file} -> {params.output_path}",
        )

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
        import pandas as pd
    except ImportError:
        return PlotQqResult(
            status="failed",
            error_message="Plotting dependencies missing. Run: pip install -e '.[plots]'",
        )

    try:
        df = pd.read_csv(params.summary_stats_file, sep="\t", low_memory=False)
        df.columns = [c.lstrip("#") for c in df.columns]

        if params.p_col not in df.columns:
            return PlotQqResult(
                status="failed",
                error_message=f"Column '{params.p_col}' not found. Available: {list(df.columns)}",
            )

        pvals = pd.to_numeric(df[params.p_col], errors="coerce").dropna()
        pvals = pvals[(pvals > 0) & (pvals <= 1)].sort_values()
        n = len(pvals)

        if n == 0:
            return PlotQqResult(status="failed", error_message="No valid p-values found")

        lambda_gc = compute_lambda_gc(pvals)

        observed = -np.log10(pvals.values)
        expected = -np.log10(np.arange(1, n + 1) / (n + 1))

        fig, ax = plt.subplots(figsize=(6, 6))
        ax.scatter(expected, observed, s=4, alpha=0.5, color="#4878CF", rasterized=True)
        ax.plot([0, expected.max()], [0, expected.max()], "r--", linewidth=1, label="Expected")

        ax.text(0.05, 0.95, f"λ_GC = {lambda_gc:.3f}", transform=ax.transAxes,
                verticalalignment="top", fontsize=10,
                bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))

        ax.set_xlabel("Expected -log10(p)")
        ax.set_ylabel("Observed -log10(p)")
        ax.set_title("QQ Plot")
        plt.tight_layout()

        Path(params.output_path).parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(params.output_path, dpi=params.dpi, bbox_inches="tight")
        plt.close()

        return PlotQqResult(
            status="success",
            output_path=params.output_path,
            lambda_gc=round(lambda_gc, 4),
            n_variants=n,
        )

    except Exception as e:
        return PlotQqResult(status="failed", error_message=str(e))
