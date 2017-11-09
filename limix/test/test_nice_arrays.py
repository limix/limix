from limix.nice_arrays import normalise_phenotype

from pandas import DataFrame
from numpy.testing import assert_allclose, assert_equal
from numpy.random import RandomState
from numpy import newaxis, stack


def test_nice_arrays_phenotype():

    random = RandomState(0)
    v = random.randn(5)
    y = [v.copy()]
    y += [v.copy()[:, newaxis]]
    y += [DataFrame(data=v.copy())]

    samples = ['sample{}'.format(i) for i in range(5)]
    data = stack([samples, v.copy()], axis=1)
    y += [DataFrame(data=data)]

    samples = ['sample{}'.format(i) for i in range(5)]
    data = stack([samples, v.copy()], axis=1)
    df = DataFrame(data=data)
    df = df.set_index(0)
    y += [df]

    for i, yi in enumerate(y):
        df = normalise_phenotype(yi, 'normal')
        assert_allclose(df.values, v[:, newaxis])
        assert_equal(df.shape[0], 5)
        assert_equal(df.shape[1], 1)
