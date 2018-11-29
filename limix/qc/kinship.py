from __future__ import division


def normalise_covariance(K, out=None):
    r"""Variance rescaling of covariance matrix ``K``.

    Let :math:`n` be the number of rows (or columns) of ``K`` and let
    :math:`m_i` be the average of the values in the i-th column.
    Gower rescaling is defined as

    .. math::

        \mathrm K \frac{n - 1}{\text{trace}(\mathrm K) - \sum m_i}.

    It works well with `Dask`_ array as log as ``out`` is ``None``.

    Notes
    -----
    The reasoning of the scaling is as follows.
    Let :math:`\mathbf g` be a vector of :math:`n` independent samples and let
    :math:`\mathrm C` be the Gower's centering matrix.
    The unbiased variance estimator is

    .. math::
        v = \sum_i \frac{(g_i-\overline g)^2}{n-1}
        =\frac{\mathrm{Tr}
        [(\mathbf g-\overline g\mathbf 1)^t(\mathbf g-\overline g\mathbf 1)]}
        {n-1}
        = \frac{\mathrm{Tr}[\mathrm C\mathbf g\mathbf g^t\mathrm C]}{n-1}

    Let :math:`\mathrm K` be the covariance matrix of :math:`\mathbf g`.
    The expectation of the unbiased variance estimator is

    .. math::

        \mathbb E[v] =
        \frac{\mathrm{Tr}[\mathrm C\mathbb E[\mathbf g\mathbf g^t]\mathrm C]}
        {n-1}
        = \frac{\mathrm{Tr}[\mathrm C\mathrm K\mathrm C]}{n-1}

    assuming that :math:`\mathbb E[g_i]=0`.
    We thus divide :math:`\mathrm K` by :math:`\mathbb E[v]` to achieve the
    desired normalisation.

    Parameters
    ----------
    K : array_like
        Covariance matrix to be normalised.
    out : array_like, optional
        Result destination. Defaults to ``None``.

    Examples
    --------
    .. doctest::

        >>> from numpy import dot, mean, zeros
        >>> from numpy.random import RandomState
        >>> from limix.qc import normalise_covariance
        >>>
        >>> random = RandomState(0)
        >>> X = random.randn(10, 10)
        >>> K = dot(X, X.T)
        >>> Z = random.multivariate_normal(zeros(10), K, 500)
        >>> print("%.3f" % mean(Z.var(1, ddof=1)))
        9.824
        >>> Kn = normalise_covariance(K)
        >>> Zn = random.multivariate_normal(zeros(10), Kn, 500)
        >>> print("%.3f" % mean(Zn.var(1, ddof=1)))
        1.008

    .. _Dask: https://dask.pydata.org/
    """
    from numpy import copyto, asarray

    K = asarray(K, float)
    c = (K.shape[0] - 1) / (K.trace() - K.mean(0).sum())
    if out is None:
        return c * K

    copyto(out, K)
    out *= c
