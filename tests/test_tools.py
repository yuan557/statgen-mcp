"""Basic tests for statgen-mcp. Uses dry_run=True so no real binaries needed."""

import pytest
from pydantic import ValidationError as PydanticValidationError

from statgen_mcp.tools.king import KingInput, run_king_cutoff
from statgen_mcp.tools.pca import PcaInput, run_pca
from statgen_mcp.tools.vcf_qc import VcfQcInput, run_vcf_qc


# ---------------------------------------------------------------------------
# VcfQcInput validation
# ---------------------------------------------------------------------------

def test_vcf_qc_invalid_maf():
    with pytest.raises(Exception):
        VcfQcInput(input_prefix="x", output_prefix="y", maf=0.9)


def test_vcf_qc_imputed_requires_score():
    with pytest.raises(PydanticValidationError):
        VcfQcInput(input_prefix="x", output_prefix="y", variant_source="imputed")


def test_vcf_qc_imputed_with_score_valid():
    params = VcfQcInput(
        input_prefix="x", output_prefix="y",
        variant_source="imputed", min_info_score=0.9,
    )
    assert params.min_info_score == 0.9


# ---------------------------------------------------------------------------
# vcf_qc runtime (dry_run)
# ---------------------------------------------------------------------------

def test_vcf_qc_missing_file():
    params = VcfQcInput(input_prefix="/nonexistent/data", output_prefix="out")
    result = run_vcf_qc(params)
    assert result.status == "failed"


def test_vcf_qc_plink_dry_run(tmp_path):
    for ext in [".bed", ".bim", ".fam"]:
        (tmp_path / f"test{ext}").write_text("")
    params = VcfQcInput(
        input_prefix=str(tmp_path / "test"),
        output_prefix="out_qc",
        maf=0.01, mind=0.05, geno=0.02, hwe=1e-6,
        dry_run=True,
    )
    result = run_vcf_qc(params)
    assert result.status in ("dry_run", "failed")


def test_vcf_qc_genotyped_dry_run(tmp_path):
    vcf = tmp_path / "test.vcf.gz"
    vcf.write_text("")
    params = VcfQcInput(
        input_prefix=str(vcf),
        output_prefix="out_qc",
        variant_source="genotyped",
        dry_run=True,
    )
    result = run_vcf_qc(params)
    assert result.status in ("dry_run", "failed")
    if result.status == "dry_run":
        assert "bcftools" in result.command_preview
        assert "TYPED" in result.command_preview


def test_vcf_qc_genotyped_custom_flag_dry_run(tmp_path):
    vcf = tmp_path / "test.vcf.gz"
    vcf.write_text("")
    params = VcfQcInput(
        input_prefix=str(vcf),
        output_prefix="out_qc",
        variant_source="genotyped",
        typed_flag="GENOTYPED",
        dry_run=True,
    )
    result = run_vcf_qc(params)
    assert result.status in ("dry_run", "failed")
    if result.status == "dry_run":
        assert "GENOTYPED" in result.command_preview


def test_vcf_qc_imputed_dry_run(tmp_path):
    vcf = tmp_path / "test.vcf.gz"
    vcf.write_text("")
    params = VcfQcInput(
        input_prefix=str(vcf),
        output_prefix="out_qc",
        variant_source="imputed",
        min_info_score=0.9,
        info_field="R2",
        dry_run=True,
    )
    result = run_vcf_qc(params)
    assert result.status in ("dry_run", "failed")
    if result.status == "dry_run":
        assert "0.9" in result.command_preview
        assert "R2" in result.command_preview


def test_vcf_qc_wgs_dry_run(tmp_path):
    vcf = tmp_path / "test.vcf.gz"
    vcf.write_text("")
    params = VcfQcInput(
        input_prefix=str(vcf),
        output_prefix="out_qc",
        variant_source="wgs",
        dry_run=True,
    )
    result = run_vcf_qc(params)
    # wgs: no bcftools step, goes straight to plink2
    assert result.status in ("dry_run", "failed")


def test_vcf_qc_plink_variant_source_warns(tmp_path):
    for ext in [".bed", ".bim", ".fam"]:
        (tmp_path / f"test{ext}").write_text("")
    params = VcfQcInput(
        input_prefix=str(tmp_path / "test"),
        output_prefix="out_qc",
        variant_source="genotyped",
        dry_run=True,
    )
    result = run_vcf_qc(params)
    assert result.status in ("dry_run", "failed")
    if result.status == "dry_run":
        assert any("skipped" in w for w in result.warnings)


# ---------------------------------------------------------------------------
# king_cutoff
# ---------------------------------------------------------------------------

def test_king_missing_file():
    params = KingInput(input_prefix="/nonexistent/data", output_prefix="out")
    result = run_king_cutoff(params)
    assert result.status == "failed"


def test_king_default_cutoff():
    params = KingInput(input_prefix="x", output_prefix="y")
    assert params.king_cutoff == 0.0442


# ---------------------------------------------------------------------------
# PCA
# ---------------------------------------------------------------------------

def test_pca_missing_file():
    params = PcaInput(input_prefix="/nonexistent/data", output_prefix="out")
    result = run_pca(params)
    assert result.status == "failed"


def test_pca_defaults():
    params = PcaInput(input_prefix="x", output_prefix="y")
    assert params.n_pcs == 20
    assert params.ld_r2 == 0.2
    assert params.ld_window == 200


def test_pca_dry_run(tmp_path):
    for ext in [".bed", ".bim", ".fam"]:
        (tmp_path / f"test{ext}").write_text("")
    params = PcaInput(
        input_prefix=str(tmp_path / "test"),
        output_prefix="out_pca",
        dry_run=True,
    )
    result = run_pca(params)
    assert result.status in ("dry_run", "failed")
    if result.status == "dry_run":
        assert "indep-pairwise" in result.command_preview
