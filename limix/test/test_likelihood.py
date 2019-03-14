import pytest
from limix._data import assert_likelihood


def test_likelihood_names():
    with pytest.raises(ValueError):
        assert_likelihood("Expon")

    valid_names = set(["normal", "bernoulli", "probit", "binomial", "poisson"])

    for n in valid_names:
        assert_likelihood(n)
