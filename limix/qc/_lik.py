def normalise_extreme_values(y, likelihood):
    from numpy import isfinite, all

    if not all(isfinite(y)):
        msg = "There are non-finite values in the the provided phenotype."
        raise ValueError(msg)

    if likelihood == "poisson":
        _poisson_normalise_extreme_values(y)
    elif likelihood == "binomial":
        _binomial_normalise_extreme_values(y)

    return y


def _poisson_normalise_extreme_values(y):
    import warnings
    from numpy import clip

    max_val = 25000.0
    if y.values.max() > max_val:
        msg = "Output values of Poisson likelihood greater"
        msg += " than {} is set to {} before applying GLMM."
        msg = msg.format(max_val, max_val)
        warnings.warn(msg)
    y.values[:] = clip(y.values, 0.0, max_val)


def _binomial_normalise_extreme_values(y):
    from numpy import minimum

    max_val = 300
    v = y.values
    ratio = v[:, 0] / v[:, 1]
    v[:, 1] = minimum(v[:, 1], max_val)
    v[:, 0] = ratio * v[:, 1]
    v[:, 0] = v[:, 0].round()
    y.values[:] = v
