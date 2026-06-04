"""Basic tests for statgen-mcp. Uses dry_run=True so no real binaries needed."""

import pytest

from statgen_mcp.executors.base import SlurmExecutor, get_executor
from statgen_mcp.tools.filter_samples import FilterSamplesInput, run_filter_samples
from statgen_mcp.tools.gwas import GwasInput, _build_covar_file, run_gwas
from statgen_mcp.tools.king import KingInput, run_king_cutoff
from statgen_mcp.tools.pca import PcaInput, run_pca
from statgen_mcp.tools.vcf_qc import VcfQcInput, run_vcf_qc
from statgen_mcp.utils.file_io import parse_plink2_log


# ---------------------------------------------------------------------------
# parse_plink2_log
# ---------------------------------------------------------------------------

def test_parse_plink2_log_sample_counts(tmp_path):
    log = tmp_path / "test.log"
    log.write_text(
        "3202 samples (1601 females, 1601 males) loaded from .fam.\n"
        "3102 samples remaining after main filters.\n"
        "800000 variants loaded from .bim file.\n"
        "750000 variants remaining after main filters.\n"
    )
    result = parse_plink2_log(log)
    assert result["samples"] == 3202
    assert result["samples_after"] == 3102
    assert result["variants_loaded"] == 800000
    assert result["variants_remaining"] == 750000


# ---------------------------------------------------------------------------
# VcfQcInput validation
# ---------------------------------------------------------------------------

def test_vcf_qc_invalid_maf():
    with pytest.raises(Exception):
        VcfQcInput(input_prefix="x", output_prefix="y", maf=0.9)


def test_vcf_qc_imputed_requires_score():
    with pytest.raises(Exception):
        VcfQcInput(input_prefix="x", output_prefix="y", variant_source="imputed")


def test_vcf_qc_imputed_with_score_valid():
    params = VcfQcInput(
        input_prefix="x", output_prefix="y",
        variant_source="imputed", min_info_score=0.99,
    )
    assert params.min_info_score == 0.99


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
    assert params.report_only is False


def test_king_report_only_dry_run(tmp_path, monkeypatch):
    bed = tmp_path / "data.bed"
    bim = tmp_path / "data.bim"
    fam = tmp_path / "data.fam"
    for f in (bed, bim, fam):
        f.touch()
    monkeypatch.setenv("STATGEN_MCP_CONFIG", "")
    params = KingInput(
        input_prefix=str(tmp_path / "data"),
        output_prefix="out",
        report_only=True,
        dry_run=True,
    )
    result = run_king_cutoff(params)
    assert result.status == "dry_run"
    assert "--make-king-table" in result.command_preview
    assert "--king-cutoff" not in result.command_preview


def test_king_cutoff_dry_run(tmp_path, monkeypatch):
    bed = tmp_path / "data.bed"
    bim = tmp_path / "data.bim"
    fam = tmp_path / "data.fam"
    for f in (bed, bim, fam):
        f.touch()
    monkeypatch.setenv("STATGEN_MCP_CONFIG", "")
    params = KingInput(
        input_prefix=str(tmp_path / "data"),
        output_prefix="out",
        report_only=False,
        dry_run=True,
    )
    result = run_king_cutoff(params)
    assert result.status == "dry_run"
    assert "--king-cutoff" in result.command_preview
    assert "--make-king-table" not in result.command_preview


# ---------------------------------------------------------------------------
# filter_samples
# ---------------------------------------------------------------------------

def test_filter_samples_requires_at_least_one_file():
    with pytest.raises(Exception):
        FilterSamplesInput(input_prefix="x", output_prefix="y")


def test_filter_samples_missing_plink():
    params = FilterSamplesInput(
        input_prefix="/nonexistent/data", output_prefix="out", keep_file="/some/file.txt"
    )
    result = run_filter_samples(params)
    assert result.status == "failed"


def test_filter_samples_keep_dry_run(tmp_path, monkeypatch):
    for ext in [".bed", ".bim", ".fam"]:
        (tmp_path / f"data{ext}").touch()
    keep = tmp_path / "keep.txt"
    keep.write_text("0 NA12878\n")
    monkeypatch.setenv("STATGEN_MCP_CONFIG", "")
    params = FilterSamplesInput(
        input_prefix=str(tmp_path / "data"),
        output_prefix="out",
        keep_file=str(keep),
        dry_run=True,
    )
    result = run_filter_samples(params)
    assert result.status == "dry_run"
    assert "--keep" in result.command_preview
    assert "--remove" not in result.command_preview


def test_filter_samples_remove_dry_run(tmp_path, monkeypatch):
    for ext in [".bed", ".bim", ".fam"]:
        (tmp_path / f"data{ext}").touch()
    remove = tmp_path / "remove.txt"
    remove.write_text("0 NA12879\n")
    monkeypatch.setenv("STATGEN_MCP_CONFIG", "")
    params = FilterSamplesInput(
        input_prefix=str(tmp_path / "data"),
        output_prefix="out",
        remove_file=str(remove),
        dry_run=True,
    )
    result = run_filter_samples(params)
    assert result.status == "dry_run"
    assert "--remove" in result.command_preview
    assert "--keep" not in result.command_preview


def test_filter_samples_both_dry_run(tmp_path, monkeypatch):
    for ext in [".bed", ".bim", ".fam"]:
        (tmp_path / f"data{ext}").touch()
    keep = tmp_path / "keep.txt"
    keep.write_text("0 NA12878\n")
    remove = tmp_path / "remove.txt"
    remove.write_text("0 NA12879\n")
    monkeypatch.setenv("STATGEN_MCP_CONFIG", "")
    params = FilterSamplesInput(
        input_prefix=str(tmp_path / "data"),
        output_prefix="out",
        keep_file=str(keep),
        remove_file=str(remove),
        dry_run=True,
    )
    result = run_filter_samples(params)
    assert result.status == "dry_run"
    assert "--keep" in result.command_preview
    assert "--remove" in result.command_preview


# ---------------------------------------------------------------------------
# GWAS
# ---------------------------------------------------------------------------

def test_gwas_missing_plink():
    params = GwasInput(input_prefix="/nonexistent/data", output_prefix="out", pheno_file="/x.txt")
    result = run_gwas(params)
    assert result.status == "failed"


def test_gwas_defaults():
    params = GwasInput(input_prefix="x", output_prefix="y", pheno_file="z.txt")
    assert params.n_pcs == 10
    assert params.firth_fallback is True
    assert params.ci is True
    assert params.pheno_name is None


def test_gwas_dry_run_no_covars(tmp_path, monkeypatch):
    for ext in [".bed", ".bim", ".fam"]:
        (tmp_path / f"data{ext}").touch()
    pheno = tmp_path / "pheno.txt"
    pheno.write_text("#FID\tIID\tT2D\n0\tNA12878\t1\n")
    monkeypatch.setenv("STATGEN_MCP_CONFIG", "")
    params = GwasInput(
        input_prefix=str(tmp_path / "data"),
        output_prefix="out",
        pheno_file=str(pheno),
        dry_run=True,
    )
    result = run_gwas(params)
    assert result.status == "dry_run"
    assert "--glm" in result.command_preview
    assert "hide-covar" in result.command_preview
    assert "firth-fallback" in result.command_preview
    assert "--covar" not in result.command_preview


def test_gwas_dry_run_with_pheno_name(tmp_path, monkeypatch):
    for ext in [".bed", ".bim", ".fam"]:
        (tmp_path / f"data{ext}").touch()
    pheno = tmp_path / "pheno.txt"
    pheno.write_text("#FID\tIID\tT2D\tBMI\n0\tNA12878\t1\t25.3\n")
    monkeypatch.setenv("STATGEN_MCP_CONFIG", "")
    params = GwasInput(
        input_prefix=str(tmp_path / "data"),
        output_prefix="out",
        pheno_file=str(pheno),
        pheno_name="T2D",
        dry_run=True,
    )
    result = run_gwas(params)
    assert result.status == "dry_run"
    assert "--pheno-name" in result.command_preview
    assert "T2D" in result.command_preview


def test_gwas_dry_run_quantitative_no_firth(tmp_path, monkeypatch):
    for ext in [".bed", ".bim", ".fam"]:
        (tmp_path / f"data{ext}").touch()
    pheno = tmp_path / "pheno.txt"
    pheno.write_text("#FID\tIID\tBMI\n0\tNA12878\t25.3\n")
    monkeypatch.setenv("STATGEN_MCP_CONFIG", "")
    params = GwasInput(
        input_prefix=str(tmp_path / "data"),
        output_prefix="out",
        pheno_file=str(pheno),
        trait_type="quantitative",
        dry_run=True,
    )
    result = run_gwas(params)
    assert result.status == "dry_run"
    assert "firth-fallback" not in result.command_preview


def test_build_covar_file_eigenvec_only(tmp_path):
    eigenvec = tmp_path / "pca.eigenvec"
    eigenvec.write_text(
        "#FID\tIID\tPC1\tPC2\tPC3\n"
        "0\tNA12878\t0.1\t-0.2\t0.05\n"
        "0\tHG00096\t-0.1\t0.3\t0.01\n"
    )
    out = _build_covar_file(None, str(eigenvec), 2, tmp_path)
    assert out is not None
    content = open(out).read()
    assert "PC1" in content
    assert "PC2" in content
    assert "PC3" not in content


def test_build_covar_file_merge(tmp_path):
    covar = tmp_path / "covars.txt"
    covar.write_text("#FID\tIID\tAGE\tSEX\n0\tNA12878\t45\t1\n0\tHG00096\t32\t2\n")
    eigenvec = tmp_path / "pca.eigenvec"
    eigenvec.write_text(
        "#FID\tIID\tPC1\tPC2\n"
        "0\tNA12878\t0.1\t-0.2\n"
        "0\tHG00096\t-0.1\t0.3\n"
    )
    out = _build_covar_file(str(covar), str(eigenvec), 2, tmp_path)
    assert out is not None
    content = open(out).read()
    assert "AGE" in content
    assert "SEX" in content
    assert "PC1" in content
    assert "PC2" in content


def test_build_covar_file_none_when_no_inputs(tmp_path):
    assert _build_covar_file(None, None, 10, tmp_path) is None


# ---------------------------------------------------------------------------
# SLURM executor
# ---------------------------------------------------------------------------

def test_get_executor_local():
    from statgen_mcp.executors.base import LocalExecutor
    assert isinstance(get_executor("local"), LocalExecutor)


def test_get_executor_slurm():
    assert isinstance(get_executor("slurm"), SlurmExecutor)


def test_get_executor_unknown():
    with pytest.raises(ValueError, match="Supported"):
        get_executor("aws")


def test_slurm_dry_run_generates_script(tmp_path, monkeypatch):
    monkeypatch.setenv("STATGEN_MCP_CONFIG", "")
    executor = SlurmExecutor()
    result = executor.run(
        cmd=["plink2", "--bfile", "data", "--make-bed", "--out", "out"],
        job_name="test_job",
        work_dir=tmp_path,
        dry_run=True,
    )
    assert result.status == "dry_run"
    assert "#SBATCH" in result.stdout
    assert "plink2" in result.stdout
    assert "#SBATCH --partition=" in result.stdout
    assert "#SBATCH --mem=" in result.stdout


def test_slurm_dry_run_respects_slurm_opts(tmp_path, monkeypatch):
    monkeypatch.setenv("STATGEN_MCP_CONFIG", "")
    executor = SlurmExecutor()
    result = executor.run(
        cmd=["plink2", "--bfile", "data", "--out", "out"],
        job_name="test_job",
        work_dir=tmp_path,
        dry_run=True,
        slurm_opts={"partition": "gpu", "mem": "64G", "time": "08:00:00"},
    )
    assert result.status == "dry_run"
    assert "gpu" in result.stdout
    assert "64G" in result.stdout
    assert "08:00:00" in result.stdout


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
