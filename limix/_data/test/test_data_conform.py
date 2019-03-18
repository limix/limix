from numpy.random import RandomState
from numpy.testing import assert_, assert_array_equal
from xarray import DataArray

from limix._data import conform_dataset


def test_data_conform_dataset():
    random = RandomState(0)

    y = random.randn(5)
    G = random.randn(5, 2)

    data = conform_dataset(y, G=G)
    y = data["y"]
    G = data["G"]
    assert_(isinstance(y, DataArray))
    assert_(isinstance(G, DataArray))
    assert_array_equal(y.shape, (5, 1))
    assert_array_equal(G.shape, (5, 2))
