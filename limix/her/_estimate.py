from .._data import conform_dataset, normalize_likelihood
from .._display import session_block, session_line
from ..qtl._assert import assert_finite


def estimate(y, lik, K, M=None, verbose=True):
    """
    Estimate the so-called narrow-sense heritability.

    It supports Normal, Bernoulli, Probit, Binomial, and Poisson phenotypes.

    Parameters
    ----------
    y : array_like
        Array of trait values of n individuals.
    lik : tuple, "normal", "bernoulli", "probit", "binomial", "poisson"
        Sample likelihood describing the residual distribution.
        Either a tuple or a string specifying the likelihood is required. The Normal,
        Bernoulli, Probit, and Poisson likelihoods can be selected by providing a
        string. Binomial likelihood on the other hand requires a tuple because of the
        number of trials: ``("binomial", array_like)``. Defaults to ``"normal"``.
    K : n×n array_like
        Sample covariance, often the so-called kinship matrix. It might be, for example,
        the estimated kinship relationship between the individuals. The provided matrix
        will be normalised as ``K = K / K.diagonal().mean()``.
    M : n×c array_like, optional
        Covariates matrix. If an array is passed, it will used as is; no normalisation
        will be performed. If ``None`` is passed, an offset will be used as the only
        covariate. Defaults to ``None``.
    verbose : bool, optional
        ``True`` to display progress and summary; ``False`` otherwise.

    Returns
    -------
    float
        Estimated heritability.

    Examples
    --------
    .. doctest::

        >>> from numpy import dot, exp, sqrt
        >>> from numpy.random import RandomState
        >>> from limix.her import estimate
        >>>
        >>> random = RandomState(0)
        >>>
        >>> G = random.randn(150, 200) / sqrt(200)
        >>> K = dot(G, G.T)
        >>> z = dot(G, random.randn(200)) + random.randn(150)
        >>> y = random.poisson(exp(z))
        >>>
        >>> print(estimate(y, 'poisson', K, verbose=False))  # doctest: +FLOAT_CMP
        0.18311439918863426

    Notes
    -----
    It will raise a ``ValueError`` exception if non-finite values are passed. Please,
    refer to the :func:`limix.qc.mean_impute` function for missing value imputation.
    """
    from numpy_sugar.linalg import economic_qs
    from numpy import pi, var, diag
    from glimix_core.glmm import GLMMExpFam
    from glimix_core.lmm import LMM

    lik = normalize_likelihood(lik)
    lik_name = lik[0]

    with session_block("Heritability analysis", disable=not verbose):

        with session_line("Normalising input...", disable=not verbose):
            data = conform_dataset(y, M=M, K=K)

        y = data["y"]
        M = data["M"]
        K = data["K"]

        assert_finite(y, M, K)

        if K is not None:
            K = K / diag(K).mean()
            QS = economic_qs(K)
        else:
            QS = None

        if lik_name == "normal":
            method = LMM(y.values, M.values, QS, restricted=True)
            method.fit(verbose=verbose)
        else:
            method = GLMMExpFam(y, lik, M.values, QS, n_int=500)
            method.fit(verbose=verbose, factr=1e6, pgtol=1e-3)

        g = method.scale * (1 - method.delta)
        e = method.scale * method.delta
        if lik_name == "bernoulli":
            e += pi * pi / 3

        v = var(method.mean())

        return g / (v + g + e)
