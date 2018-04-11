from numpy import clip, minimum, isfinite, all


def _poisson_normalise_extreme_values(y):
    max_val = 25000.
    y.values[:] = clip(y.values, 0., max_val)


def _binomial_normalise_extreme_values(y):
    max_val = 300
    v = y.values
    ratio = v[:, 0] / v[:, 1]
    v[:, 1] = minimum(v[:, 1], max_val)
    v[:, 0] = ratio * v[:, 1]
    v[:, 0] = v[:, 0].round()
    y.values[:] = v


def normalise_extreme_values(y, likelihood):
    if not all(isfinite(y)):
        msg = "There are non-finite values in the the provided phenotype."
        raise ValueError(msg)

    if likelihood == 'poisson':
        _poisson_normalise_extreme_values(y)
    elif likelihood == 'binomial':
        _binomial_normalise_extreme_values(y)
