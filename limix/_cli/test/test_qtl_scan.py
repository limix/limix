import pytest
from click.testing import CliRunner
from numpy.testing import assert_equal

from limix import cli


def _urls(suffixes):
    base = "http://rest.s3for.me/limix/"
    return [base + s for s in suffixes]


@pytest.mark.remfiles(_urls(["trait.csv"]))
def test_cli_qtl_scan_no_variant(remfiles):
    remfiles.chdir()
    invoke = CliRunner().invoke
    r = invoke(cli, ["qtl", "scan", "--trait", "trait.csv"])
    assert_equal("Error: no variant has been specified." in r.stdout, True)
    assert_equal(r.exit_code, 2)


@pytest.mark.remfiles(_urls(["trait.csv", "plink.bed"]))
def test_cli_qtl_scan_trait_bed(remfiles):
    remfiles.chdir()
    invoke = CliRunner().invoke
    r = invoke(cli, ["qtl", "scan", "--bfile", "plink", "--trait", "trait.csv"])
    assert_equal(type(r.exception), FileNotFoundError)
    assert_equal(r.exit_code, 1)


@pytest.mark.remfiles(_urls(["trait.csv", "plink.bed", "plink.bim"]))
def test_cli_qtl_scan_trait_bed_bim(remfiles):
    remfiles.chdir()
    invoke = CliRunner().invoke
    r = invoke(cli, ["qtl", "scan", "--bfile", "plink", "--trait", "trait.csv"])
    assert_equal(type(r.exception), FileNotFoundError)
    assert_equal(r.exit_code, 1)


@pytest.mark.remfiles(_urls(["trait.csv", "plink.bed", "plink.bim", "plink.fam"]))
def test_cli_qtl_scan_trait_bed_bim_fam(remfiles):
    remfiles.chdir()
    runner = CliRunner()
    r = runner.invoke(cli, ["qtl", "scan", "--bfile", "plink", "--trait", "trait.csv"])
    assert_equal(r.exit_code, 0)


@pytest.mark.remfiles(
    _urls(
        [
            "rel/plink2.rel.bin",
            "rel/plink2.rel.id",
            "grm-bin/plink.grm.N.bin",
            "grm-bin/plink.grm.bin",
            "grm-bin/plink.grm.id",
            "trait.csv",
        ]
    )
)
def test_cli_qtl_scan_grm_rel(remfiles):
    remfiles.chdir()
    runner = CliRunner()
    r = runner.invoke(
        cli,
        [
            "qtl",
            "scan",
            "--bfile",
            "plink",
            "--rel",
            "plink2.rel.bin",
            "--grm",
            "plink.grm.bin",
            "--trait",
            "trait.csv",
        ],
    )
    msg = "Error: The options [--grm, --rel] are mutually exclusive."
    assert_equal(msg in r.stdout, True)
    assert_equal(r.exit_code, 1)


@pytest.mark.remfiles(_urls(["plink.bed", "plink.bim", "plink.fam", "trait.csv"]))
def test_cli_qtl_scan_traits_not_found(remfiles):
    remfiles.chdir()
    runner = CliRunner()
    r = runner.invoke(
        cli,
        [
            "qtl",
            "scan",
            "--bfile",
            "plink",
            "--trait",
            "trait.csv",
            "--method=st",
            "--trait-name",
            "pinto,pepeca",
        ],
    )
    ok = "Error: not all specified traits have been found." in r.stdout
    assert_equal(ok, True)
    assert_equal(r.exit_code, 2)
