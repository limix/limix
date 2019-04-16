import limix
from numpy.testing import assert_equal


def test_threads():
    limix.threads.set_max_nthreads(5)
    assert_equal(limix.threads.get_max_nthreads(), 5)
