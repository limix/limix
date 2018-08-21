from os.path import dirname, join, realpath

from numpy.testing import assert_allclose, assert_equal

from limix import io


def test_io_hdf5():
    dir_path = dirname(realpath(__file__))
    with io.hdf5.fetcher(join(dir_path, "data", "hdf5", "data.h5")) as df:
        X = df.fetch("/genotype/chrom20/SNP")
        assert_allclose(X[0, 0], [1.0])


def test_io_csv():
    dir_path = dirname(realpath(__file__))
    fn = join(dir_path, "data", "csv", "pheno.csv")
    assert_equal(io.csv.read(fn)["attr1"][0], "string")


def test_io_gen():
    dir_path = dirname(realpath(__file__))
    df = io.gen.read(join(dir_path, "data", "gen", "example"))
    assert_equal(df["sample"]["sample_id"][1], "1A1")
    assert_equal(df["sample"]["age"][0], 4)
    assert_allclose(df["genotype"]["1A4"]["AB"][0], 0.0207)
