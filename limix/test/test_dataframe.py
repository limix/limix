from limix.dataframe import infer_samples_index
from limix.dataframe import assert_compatible_samples_arrays
import pytest
from numpy import array
from numpy.random import RandomState
from numpy.testing import assert_, assert_array_equal
from pandas import DataFrame, Index


def test_dataframe_arrays_compatibility():
    y = array([-1.2, 3.4, 0.1])
    assert_(assert_compatible_samples_arrays([]))
    assert_(assert_compatible_samples_arrays(y))
    assert_(assert_compatible_samples_arrays([y]))

    random = RandomState(0)

    K = random.randn(3, 4)
    K = K.dot(K.T)
    assert_(assert_compatible_samples_arrays([y, K]))

    M = random.randn(3, 2)
    assert_(assert_compatible_samples_arrays([y, K, M]))
    Me = random.randn(2, 2)
    assert_(not assert_compatible_samples_arrays([y, K, Me]))


def test_dataframe_samples_index_infer():
    y = array([-1.2, 3.4, 0.1])

    random = RandomState(0)

    K = random.randn(3, 4)
    K = K.dot(K.T)

    M = random.randn(3, 2)
    Me = random.randn(2, 2)

    with pytest.raises(ValueError):
        infer_samples_index([y, K, Me])

    candidates = Index(['candidate0', 'candidate1', 'candidate2'])
    assert_array_equal(infer_samples_index(y), candidates)
    assert_array_equal(infer_samples_index([y]), candidates)
    assert_array_equal(infer_samples_index([y, K]), candidates)
    assert_array_equal(infer_samples_index([y, K, M]), candidates)

    df_y = DataFrame(y, index=['c0', 'c1', 'c2'])
    df_K = DataFrame(K, index=['c0', 'c1', 'c2'])
    df_Me = DataFrame(Me, index=['c0', 'c1'])
    c = Index(['c0', 'c1'])
    assert_array_equal(infer_samples_index([df_y, df_K, df_Me]), c)

    cindex = Index(['c0', 'c1', 'c2'])
    assert_array_equal(infer_samples_index([df_y, df_K]), cindex)
    assert_array_equal(infer_samples_index([df_y, df_K, M]), cindex)

    df_K = DataFrame(K, index=['c2', 'c0', 'c1'])
    with pytest.raises(ValueError):
        infer_samples_index([df_y, df_K, M])
