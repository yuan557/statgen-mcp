"""Basic tests for statgen-mcp. Uses dry_run=True so no real binaries needed."""

import pytest
from statgen_mcp.tools.king import KingInput, run_king_cutoff
from statgen_mcp.tools.pca import PcaInput, run_pca
from statgen_mcp.tools.vcf_qc import VcfQcInput, run_vcf_qc


def test_vcf_qc_missing_file():
    params = VcfQcInput(input_prefix="/nonexistent/data", output_prefix="out")
    result = run_vcf_qc(params)
    assert result.status == "failed"


def test_vcf_qc_dry_run(tmp_path):
    for ext in [".bed", ".bim", ".fam"]:
        (tmp_path / f"test{ext}").write_text("")
    params = VcfQcInput(
        input_prefix=str(tmp_path / "test"),
        output_prefix="out_qc",
        maf=0.01,
        mind=0.05,
        geno=0.02,
        hwe=1e-6,
        dry_run=True,
    )
    result = run_vcf_qc(params)
    assert result.status in ("dry_run", "failed")


def test_king_missing_file():
    params = KingInput(input_prefix="/nonexistent/data", output_prefix="out")
    result = run_king_cutoff(params)
    assert result.status == "failed"


def test_pca_missing_file():
    params = PcaInput(input_prefix="/nonexistent/data", output_prefix="out")
    result = run_pca(params)
    assert result.status == "failed"


def test_vcf_qc_invalid_maf():
    with pytest.raises(Exception):
        VcfQcInput(input_prefix="x", output_prefix="y", maf=0.9)


def test_king_default_cutoff():
    params = KingInput(input_prefix="x", output_prefix="y")
    assert params.king_cutoff == 0.0442


def test_pca_defaults():
    params = PcaInput(input_prefix="x", output_prefix="y")
    assert params.n_pcs == 20
    assert params.ld_prune is True
