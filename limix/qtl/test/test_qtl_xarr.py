import limix
import pytest
from numpy.random import RandomState
import xarray as xr
from numpy.testing import assert_allclose


def test_qtl_xarr():
    with limix.example.file_example("xarr.hdf5.bz2") as filepath:
        filepath = limix.sh.extract(filepath, verbose=False)
        sample_ids = limix.io.hdf5.fetch(filepath, "/foo/chr1/col_header/sample_ids")

        rsid = dict()
        for i in range(1, 3):
            rsid[i] = limix.io.hdf5.fetch(
                filepath, "/foo/chr{}/row_header/rsid".format(i)
            )

        G = []
        for i in range(1, 3):
            g = xr.open_dataset(filepath, "/foo/chr{}".format(i))["matrix"]
            g = g.rename({g.dims[0]: "snps", g.dims[1]: "samples"})
            g = g.T
            g["samples"] = sample_ids
            g["snps"] = rsid[i]
            g["chr"] = ("snps", [i] * g.shape[1])
            G.append(g)

        G = xr.concat(G, dim="snps")

        K = limix.stats.linear_kinship(G)

    k = [
        [4.316150754626438, -0.12182214897716158],
        [-0.12182214897716158, 0.25268948339191105],
    ]
    assert_allclose(K[:2][:, :2], k)
    random = RandomState(0)
    y = random.randn(10)

    with pytest.raises(ValueError):
        limix.qtl.scan(G, y, "normal", K)

    G = G.rename({"samples": "sample"})
    limix.qtl.scan(G, y, "normal", K)
