from click import UsageError
from click.testing import CliRunner
from numpy.testing import assert_equal

from limix import cli
from limix.sh import download


def test_cli_qtl_scan():

    runner = CliRunner()
    with runner.isolated_filesystem():
        invoke = runner.invoke
        download("http://rest.s3for.me/limix/example.pphe", verbose=False)
        result = invoke(cli, ["qtl", "scan", "--pheno", "example.pphe"])
        assert_equal(result.exit_code, 2)

        download("http://rest.s3for.me/limix/plink.bed", verbose=False)
        result = invoke(
            cli, ["qtl", "scan", "--bfile", "plink", "--pheno", "example.pphe"]
        )
        assert_equal(type(result.exception), FileNotFoundError)
        assert_equal(result.exit_code, 1)

        download("http://rest.s3for.me/limix/plink.bim", verbose=False)
        result = invoke(
            cli, ["qtl", "scan", "--bfile", "plink", "--pheno", "example.pphe"]
        )
        assert_equal(type(result.exception), FileNotFoundError)
        assert_equal(result.exit_code, 1)

        download("http://rest.s3for.me/limix/plink.fam", verbose=False)
        result = invoke(
            cli, ["qtl", "scan", "--bfile", "plink", "--pheno", "example.pphe"]
        )
        assert_equal(result.exit_code, 0)

        result = invoke(cli, ["qtl", "scan", "--pheno", "example.pphe"])
        assert_equal("Error: No variant has been specified." in result.stdout, True)
        assert_equal(result.exit_code, 2)

        result = invoke(cli, ["qtl", "scan", "--bfile", "plink"])
        assert_equal("Error: No phenotype has been specified." in result.stdout, True)
        assert_equal(result.exit_code, 2)
