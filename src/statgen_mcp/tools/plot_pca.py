"""PCA scatter plot coloured by population label."""

from __future__ import annotations

from pathlib import Path

from pydantic import BaseModel, Field

from statgen_mcp.utils.validation import ValidationError, require_file


DEFAULT_COLORS = {
    "EUR": "#4878CF",
    "AFR": "#D65F5F",
    "EAS": "#6ACC65",
    "AMR": "#B47CC7",
    "SAS": "#C4AD66",
}


class PlotPcaInput(BaseModel):
    eigenvec_file: str = Field(description="PLINK2 .eigenvec file from run_pca")
    output_path: str = Field(description="Output PNG file path")
    label_file: str | None = Field(
        default=None,
        description="File with sample labels for colouring (must contain IID and label_col columns)",
    )
    label_col: str = Field(
        default="Superpopulation",
        description="Column name in label_file to colour points by",
    )
    label_id_col: str = Field(
        default="SampleID",
        description="Column name in label_file that matches IID in eigenvec",
    )
    pc_pairs: list[list[int]] = Field(
        default=[[1, 2], [3, 4]],
        description="PC pairs to plot, e.g. [[1,2],[3,4]] produces two panels",
    )
    point_size: float = Field(default=8.0, description="Scatter point size")
    alpha: float = Field(default=0.7, description="Scatter point alpha")
    dpi: int = Field(default=150, description="Output image DPI")
    dry_run: bool = Field(default=False)


class PlotPcaResult(BaseModel):
    status: str
    output_path: str = ""
    n_samples: int | None = None
    n_labels: int | None = None
    error_message: str | None = None
    command_preview: str | None = None


def run_plot_pca(params: PlotPcaInput) -> PlotPcaResult:
    try:
        require_file(params.eigenvec_file)
        if params.label_file:
            require_file(params.label_file)
    except ValidationError as e:
        return PlotPcaResult(status="failed", error_message=str(e))

    if params.dry_run:
        preview = f"plot_pca: {params.eigenvec_file}"
        if params.label_file:
            preview += f" coloured by {params.label_col} from {params.label_file}"
        preview += f" → {params.output_path}"
        return PlotPcaResult(status="dry_run", command_preview=preview)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.patches as mpatches
        import matplotlib.pyplot as plt
        import pandas as pd
    except ImportError:
        return PlotPcaResult(
            status="failed",
            error_message="Plotting dependencies missing. Run: pip install -e '.[plots]'",
        )

    try:
        pcs = pd.read_csv(params.eigenvec_file, sep="\t", comment="#", header=0)
        pcs.columns = ["FID", "IID"] + [f"PC{i}" for i in range(1, len(pcs.columns) - 1)]

        n_samples = len(pcs)
        n_labels = None

        if params.label_file:
            labels = pd.read_csv(params.label_file, sep=r"\s+")
            pcs = pcs.merge(
                labels[[params.label_id_col, params.label_col]],
                left_on="IID", right_on=params.label_id_col, how="left",
            )
            n_labels = pcs[params.label_col].nunique()

        n_panels = len(params.pc_pairs)
        fig, axes = plt.subplots(1, n_panels, figsize=(7 * n_panels, 6))
        if n_panels == 1:
            axes = [axes]

        unique_labels = (
            sorted(pcs[params.label_col].dropna().unique())
            if params.label_file
            else None
        )

        for ax, pair in zip(axes, params.pc_pairs):
            xpc, ypc = f"PC{pair[0]}", f"PC{pair[1]}"
            if params.label_file and params.label_col in pcs.columns:
                for label, grp in pcs.groupby(params.label_col):
                    color = DEFAULT_COLORS.get(str(label), None)
                    ax.scatter(
                        grp[xpc], grp[ypc],
                        c=color, s=params.point_size, alpha=params.alpha, label=label,
                    )
            else:
                ax.scatter(
                    pcs[xpc], pcs[ypc],
                    s=params.point_size, alpha=params.alpha, color="#4878CF",
                )
            ax.set_xlabel(xpc)
            ax.set_ylabel(ypc)
            ax.set_title(f"{xpc} vs {ypc}")

        if unique_labels is not None:
            handles = [
                mpatches.Patch(color=DEFAULT_COLORS.get(l, "grey"), label=l)
                for l in unique_labels
            ]
            fig.legend(handles=handles, loc="lower center", ncol=min(len(unique_labels), 5),
                       bbox_to_anchor=(0.5, -0.05))

        plt.tight_layout()
        Path(params.output_path).parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(params.output_path, dpi=params.dpi, bbox_inches="tight")
        plt.close()

        return PlotPcaResult(
            status="success",
            output_path=params.output_path,
            n_samples=n_samples,
            n_labels=n_labels,
        )

    except Exception as e:
        return PlotPcaResult(status="failed", error_message=str(e))
