import pytest
from limix.likelihood import assert_likelihood_name


def test_likelihood_names():
    with pytest.raises(ValueError):
        assert_likelihood_name("Expon")

    valid_names = set(["normal", "bernoulli", "probit", "binomial", "poisson"])

    for n in valid_names:
        assert_likelihood_name(n)
