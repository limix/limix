from numpy.testing import assert_allclose, assert_equal

from limix import file_example, io
from limix.sh import extract


def test_io_hdf5():
    with file_example("data.h5.bz2") as filepath:
        filepath = extract(filepath, verbose=False)
        with io.hdf5.fetcher(filepath) as df:
            X = df.fetch("/genotype/chrom20/SNP")
            assert_allclose(X[0, 0], [1.0])


def test_io_csv():
    with file_example("pheno.csv") as filepath:
        assert_equal(io.csv.read(filepath, verbose=False)["attr1"][0], "string")


def test_io_gen():
    with file_example(["example.gen", "example.sample"]) as filepaths:
        df = io.gen.read(filepaths[0][:-4], verbose=False)
        assert_equal(df["sample"]["sample_id"][1], "1A1")
        assert_equal(df["sample"]["age"][0], 4)
        assert_allclose(df["genotype"]["1A4"]["AB"][0], 0.0207)
