from click.testing import CliRunner
from numpy.testing import assert_equal

from limix import cli
from limix.sh import download


def test_cli_qtl_scan():

    runner = CliRunner()
    with runner.isolated_filesystem():
        invoke = runner.invoke
        download("http://rest.s3for.me/limix/trait.csv", verbose=False)
        result = invoke(cli, ["qtl", "scan", "--trait", "trait.csv"])
        assert_equal(result.exit_code, 2)

        download("http://rest.s3for.me/limix/plink.bed", verbose=False)
        result = invoke(
            cli, ["qtl", "scan", "--bfile", "plink", "--trait", "trait.csv"]
        )
        assert_equal(type(result.exception), FileNotFoundError)
        assert_equal(result.exit_code, 1)

        download("http://rest.s3for.me/limix/plink.bim", verbose=False)
        result = invoke(
            cli, ["qtl", "scan", "--bfile", "plink", "--trait", "trait.csv"]
        )
        assert_equal(type(result.exception), FileNotFoundError)
        assert_equal(result.exit_code, 1)

        download("http://rest.s3for.me/limix/plink.fam", verbose=False)
        result = invoke(
            cli, ["qtl", "scan", "--bfile", "plink", "--trait", "trait.csv"]
        )
        assert_equal(result.exit_code, 0)

        result = invoke(cli, ["qtl", "scan", "--trait", "trait.csv"])
        assert_equal("Error: no variant has been specified." in result.stdout, True)
        assert_equal(result.exit_code, 2)

        result = invoke(cli, ["qtl", "scan", "--bfile", "plink"])
        assert_equal("Error: no trait has been specified." in result.stdout, True)
        assert_equal(result.exit_code, 2)

        download("http://rest.s3for.me/limix/rel/plink2.rel.bin", verbose=False)
        download("http://rest.s3for.me/limix/rel/plink2.rel.id", verbose=False)

        download("http://rest.s3for.me/limix/grm-bin/plink.grm.N.bin", verbose=False)
        download("http://rest.s3for.me/limix/grm-bin/plink.grm.bin", verbose=False)
        download("http://rest.s3for.me/limix/grm-bin/plink.grm.id", verbose=False)
        result = invoke(
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
        assert_equal(
            "Error: The options [--grm, --rel] are mutually exclusive."
            in result.stdout,
            True,
        )
        assert_equal(result.exit_code, 1)
        result = invoke(
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
        ok = "Error: not all specified traits have been found." in result.stdout
        assert_equal(ok, True)
        assert_equal(result.exit_code, 2)
