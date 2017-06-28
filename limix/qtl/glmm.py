from __future__ import division

from time import time

from numpy import abs as npy_abs
from numpy import asarray as npy_asarray
from numpy import diag, ones, sqrt
from numpy_sugar.linalg import economic_qs
from scipy.stats import chi2

from glimix_core.glmm import GLMM
from limix.qtl.lmm import LMM
from limix.util import asarray

from .qtl_model import QTLModel


def qtl_test_glmm(G, y, lik, K, M=None, verbose=True):
    r"""Single-variant association testing via generalised linear mixed models.

    Parameters
    ----------
    G : array_like
        `N` individuals by `S` SNPs.
    y : (tuple, array_like)
        Either a tuple of two arrays of `N` individuals each (Binomial
        phenotypes) or an array of `N` individuals (Poisson or Bernoulli
        phenotypes). It does not support missing values yet.
    lik : {'bernoulli', 'binomial', 'poisson'}
        Sample likelihood describing the residual distribution.
    K : array_like
        `N` by `N` covariance matrix (e.g., kinship coefficients).
    M : (array_like, optional)
        `N` individuals by `D` covariates.
        By default, ``M`` is a (`N`, `1`) array of ones.
    verbose : (bool, optional)
        if ``True``, details such as runtime are displayed.

    Returns
    -------
        :class:`limix.qtl.LMM`: LIMIX LMM object

    Examples
    --------
    .. doctest::

        >>> from numpy import dot, exp, sqrt
        >>> from numpy.random import RandomState
        >>> from limix.qtl import qtl_test_glmm
        >>>
        >>> random = RandomState(0)
        >>>
        >>> G = random.randn(250, 500) / sqrt(500)
        >>> beta = 0.01 * random.randn(500)
        >>>
        >>> z = dot(G, beta) + 0.1 * random.randn(250)
        >>> z += dot(G[:, 0], 1) # causal SNP
        >>>
        >>> y = random.poisson(exp(z))
        >>>
        >>> candidates = G[:, :5]
        >>> K = dot(G[:, 5:], G[:, 5:].T)
        >>> model = qtl_test_glmm(candidates, y, 'poisson', K, verbose=False)
        >>>
        >>> print(model.variant_pvalues)
        [ 0.0694  0.3336  0.5899  0.7388  0.7796]
    """

    G = asarray(G)

    if M is None:
        M = ones((G.shape[0], 1))
    else:
        M = asarray(M)

    K = asarray(K)

    if isinstance(y, (tuple, list)):
        y = tuple([npy_asarray(p, float) for p in y])
    else:
        y = npy_asarray(y, float)

    start = time()
    QS = economic_qs(K)
    glmm = GLMM(y, lik, M, QS)
    glmm.feed().maximize(progress=verbose)

    # extract stuff from glmm
    eta = glmm._site.eta
    tau = glmm._site.tau
    scale = float(glmm.scale)
    delta = float(glmm.delta)

    # define useful quantities
    mu = eta / tau
    var = 1. / tau
    s2_g = scale * (1 - delta)
    tR = s2_g * K + diag(var - var.min() + 1e-4)

    start = time()
    lmm = LMM(snps=G, pheno=mu, K=tR, covs=M, verbose=verbose)

    return QTLModel_GLMM(lmm)


class QTLModel_GLMM(QTLModel):
    def __init__(self, lmm):
        self._lmm = lmm

    @property
    def null_lml(self):
        pass

    @property
    def alt_lmls(self):
        pass

    @property
    def variant_effsizes(self):
        return self._lmm.getBetaSNP().ravel()

    @property
    def variant_effsizes_se(self):
        return self._lmm.getBetaSNPste().ravel()

    @property
    def variant_pvalues(self):
        return self._lmm.getPv().ravel()

    @property
    def null_covariant_effsizes(self):
        pass

    @property
    def alt_covariant_effsizes(self):
        return dict(offset=ones(5) * 1.5)

    def __str__(self):
        from pandas import DataFrame

        data = dict(
            effsizes=self.variant_effsizes,
            effsizes_se=self.variant_effsizes_se,
            pvalues=self.variant_pvalues)

        variant_msg = str(DataFrame(data=data).describe())

        data = dict(self.alt_covariant_effsizes)

        covariant_msg = str(DataFrame(data=data).describe())

        msg = 'variants\n' + variant_msg
        msg += '\n\ncovariants\n' + covariant_msg

        return msg
