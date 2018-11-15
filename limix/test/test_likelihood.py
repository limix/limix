import pytest
from limix._data._lik import check_likelihood_name


def test_likelihood_names():
    with pytest.raises(ValueError):
        check_likelihood_name("Expon")

    valid_names = set(["normal", "bernoulli", "probit", "binomial", "poisson"])

    for n in valid_names:
        check_likelihood_name(n)
