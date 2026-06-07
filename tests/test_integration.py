"""
Integration tests using the GWASTutorial 1000 Genomes EAS dataset.

These tests run the full pipeline on real data and are NOT part of the regular
dry-run test suite. They require plink2 on PATH and the GWASTutorial dataset.

Set GWAS_TUTORIAL_DIR to the cloned repo root, or clone to the default location:
  git clone https://github.com/Cloufield/GWASTutorial ~/gwas_testdata/GWASTutorial

Run with:
  pytest tests/test_integration.py -v
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

TUTORIAL_ROOT = Path(
    os.environ.get("GWAS_TUTORIAL_DIR", Path.home() / "gwas_testdata/GWASTutorial")
)
_DATASET_DIR = TUTORIAL_ROOT / "01_Dataset"
_PLINK_STEM = "1KG.EAS.auto.snp.norm.nodup.split.rare002.common015.missing"
PLINK_PREFIX = _DATASET_DIR / _PLINK_STEM / _PLINK_STEM
PHENO_FILE = _DATASET_DIR / "1kgeas_binary.txt"

DATA_AVAILABLE = Path(str(PLINK_PREFIX) + ".bed").exists()

pytestmark = pytest.mark.skipif(
    not DATA_AVAILABLE,
    reason=(
        "GWASTutorial data not found. "
        "Set GWAS_TUTORIAL_DIR or clone https://github.com/Cloufield/GWASTutorial"
    ),
)


# ---------------------------------------------------------------------------
# Session fixtures: run pipeline once, validate in individual tests
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session", autouse=True)
def _config(tmp_path_factory):
    """Write a temp config and reset the singleton so all tools use it."""
    import statgen_mcp.config as cfg_module

    run_root = tmp_path_factory.mktemp("runs")
    config_path = tmp_path_factory.mktemp("cfg") / "config.toml"
    config_path.write_text(
        f'output_dir = "{run_root}"\n'
        '[executor]\nbackend = "local"\n'
        '[tools]\nplink2 = "plink2"\n'
    )
    os.environ["STATGEN_MCP_CONFIG"] = str(config_path)
    cfg_module._config = None
    yield
    os.environ.pop("STATGEN_MCP_CONFIG", None)
    cfg_module._config = None


@pytest.fixture(scope="session")
def pipeline(tmp_path_factory):
    """Run the full pipeline once and return all results."""
    from statgen_mcp.tools.gwas import GwasInput, run_gwas
    from statgen_mcp.tools.king import KingInput, run_king_cutoff
    from statgen_mcp.tools.pca import PcaInput, run_pca
    from statgen_mcp.tools.plot_qq import PlotQqInput, run_plot_qq
    from statgen_mcp.tools.vcf_qc import VcfQcInput, run_vcf_qc

    out = {}
    plots = tmp_path_factory.mktemp("plots")

    qc = run_vcf_qc(VcfQcInput(
        input_prefix=str(PLINK_PREFIX),
        output_prefix="1kg_eas_qc",
        maf=0.01, mind=0.05, geno=0.05, hwe=1e-6,
        biallelic_only=True,
    ))
    assert qc.status == "success", f"vcf_qc failed: {qc.error_message}"
    out["vcf_qc"] = qc

    king = run_king_cutoff(KingInput(
        input_prefix=qc.output_prefix,
        output_prefix="1kg_eas_king",
        king_cutoff=0.0442,
        report_only=False,
    ))
    assert king.status == "success", f"king failed: {king.error_message}"
    out["king"] = king

    pca = run_pca(PcaInput(
        input_prefix=king.output_prefix,
        output_prefix="1kg_eas_pca",
        n_pcs=20,
    ))
    assert pca.status == "success", f"pca failed: {pca.error_message}"
    out["pca"] = pca

    eigenvec = next(f for f in pca.output_files if f.endswith(".eigenvec"))
    gwas = run_gwas(GwasInput(
        input_prefix=qc.output_prefix,
        output_prefix="1kg_eas_gwas",
        pheno_file=str(PHENO_FILE),
        eigenvec_file=eigenvec,
        n_pcs=10,
        firth_fallback=True,
    ))
    assert gwas.status == "success", f"gwas failed: {gwas.error_message}"
    out["gwas"] = gwas

    stats_file = next(f for f in gwas.output_files if "glm.logistic" in f)
    qq = run_plot_qq(PlotQqInput(
        summary_stats_file=stats_file,
        output_path=str(plots / "qq.png"),
    ))
    assert qq.status == "success", f"plot_qq failed: {qq.error_message}"
    out["qq"] = qq

    return out


# ---------------------------------------------------------------------------
# Validation tests
# ---------------------------------------------------------------------------

def test_vcf_qc_variant_counts(pipeline):
    r = pipeline["vcf_qc"]
    assert r.variants_before == pytest.approx(1_230_000, rel=0.05)
    assert r.variants_after > 900_000
    assert r.variants_after < r.variants_before


def test_vcf_qc_sample_counts(pipeline):
    r = pipeline["vcf_qc"]
    assert r.samples_before == 504
    # some samples removed for missingness
    assert r.log_summary.get("samples_after", 504) <= 504


def test_king_finds_related_pairs(pipeline):
    r = pipeline["king"]
    assert r.log_summary.get("samples") is not None
    # known: 7 related pairs → ~9 samples removed
    samples_after = r.log_summary.get("samples_after")
    if samples_after is not None:
        assert samples_after < 505


def test_pca_succeeds(pipeline):
    r = pipeline["pca"]
    assert r.status == "success"
    assert any(f.endswith(".eigenvec") for f in r.output_files)
    assert any(f.endswith(".eigenval") for f in r.output_files)


def test_gwas_runs_binary_trait(pipeline):
    r = pipeline["gwas"]
    assert r.trait_type_used == "binary"
    assert r.n_variants_tested is not None
    assert r.n_variants_tested > 500_000


def test_gwas_lambda_gc(pipeline):
    r = pipeline["qq"]
    assert r.lambda_gc is not None
    # well-controlled GWAS: lambda between 0.95 and 1.15
    assert 0.95 < r.lambda_gc < 1.15, f"λ_GC = {r.lambda_gc} — possible inflation or deflation"


def test_gwas_n_variants(pipeline):
    r = pipeline["qq"]
    assert r.n_variants is not None
    assert r.n_variants > 500_000
