import numpy as np
import pytest
import xarray as xr
from numpy.testing import assert_array_equal as assert_eq, assert_string_equal

from limix._data import asarray


def test_asarray():
    from limix._data._conf import CONF

    random = np.random.RandomState(0)
    arr = random.randn(2, 3)

    def _assert_target(x, target):
        assert_eq(x.dims, CONF["data_dims"][target])
        assert_string_equal(x.name, target)

    _assert_target(asarray(arr, "trait", {"row": "sample"}), "trait")
    _assert_target(asarray(arr, "trait", {"col": "sample"}), "trait")
    _assert_target(asarray(arr, "trait", {"row": "trait"}), "trait")
    _assert_target(asarray(arr, "trait", {"col": "trait"}), "trait")
    _assert_target(asarray(arr, "trait", {"row": "trait", "col": "sample"}), "trait")
    _assert_target(asarray(arr, "trait", {"row": "trait", 1: "sample"}), "trait")
    _assert_target(asarray(arr, "trait", {0: "trait", 1: "sample"}), "trait")
    _assert_target(asarray(arr, "trait", ["trait"]), "trait")
    _assert_target(asarray(arr, "trait", ["sample"]), "trait")
    _assert_target(asarray(arr, "trait", ["sample", "trait"]), "trait")
    _assert_target(asarray(arr, "trait", ["sample", "wrong"]), "trait")

    with pytest.raises(ValueError):
        _assert_target(asarray(arr, "wrong"), "trait")

    with pytest.raises(ValueError):
        _assert_target(asarray(arr, "trait", ["wrong1", "wrong2"]), "trait")

    _assert_target(
        asarray(xr.DataArray(arr, dims=["sample", "trait"]), "trait"), "trait"
    )
    _assert_target(
        asarray(xr.DataArray(arr, dims=["sample", "dim2"]), "trait"), "trait"
    )
    _assert_target(asarray(xr.DataArray(arr, dims=["dim1", "trait"]), "trait"), "trait")
    _assert_target(
        asarray(
            xr.DataArray(arr, dims=["dim1", "trait"]), "trait", ["sample", "trait"]
        ),
        "trait",
    )
    _assert_target(
        asarray(xr.DataArray(arr, dims=["dim1", "trait"]), "trait", {0: "sample"}),
        "trait",
    )
    _assert_target(
        asarray(xr.DataArray(arr, dims=["dim1", "trait"]), "trait", {1: "sample"}),
        "trait",
    )

    with pytest.raises(ValueError):
        _assert_target(
            asarray(xr.DataArray(arr), "trait", {0: "sample", 1: "sample"}), "trait"
        )

    with pytest.raises(ValueError):
        _assert_target(
            asarray(xr.DataArray(arr), "trait", ["sample", "sample"]), "trait"
        )
