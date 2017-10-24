r"""
**********************
Single-variant testing
**********************


Introduction
^^^^^^^^^^^^

A **linear mixed model** (LMM) can be described as

.. math::

    \mathbf y = \mathrm M\boldsymbol\beta_0 + \mathrm X\mathbf u
    + \boldsymbol\epsilon,

where :math:`\mathbf u \sim \mathcal N(\mathbf 0, \sigma_u^2\mathrm I)` is a
vector of random effects and :math:`\epsilon_i` are iid Normal random
variables with zero-mean and variance :math:`\sigma_e^2` each.
Covariates are defined by the columns of :math:`\mathrm M`, and
:math:`\mathrm X` commonly contain all genetic variants of each sample.

The outcome-vector is thus distributed according to

.. math::

    \mathbf y \sim \mathcal N(\mathrm M\boldsymbol\beta_0,
                              \sigma_u^2 \mathrm X \mathrm X^{\intercal}
                              + \sigma_e^2\mathrm I).

The parameters :math:`\boldsymbol\beta_0`, :math:`\sigma_u`, and
:math:`\sigma_{\epsilon}` are estimated via the maximum likelihood estimation
(MLE) approach under the null hypothesis just defined.

The alternative hypothesis for single-variant testing consists in the addition
of a fixed-effect size :math:`\beta_1`:

.. math::

    \mathbf y = \mathrm M\boldsymbol\beta_1 + \mathbf g\beta_1
        + \mathrm X\mathbf u + \boldsymbol\epsilon.

The parameter :math:`\beta_1` multiplies a given vector :math:`\mathbf g`,
typically representing a genetic locus of interest.
The parameters :math:`\boldsymbol\beta_0`, :math:`\beta_1`,
:math:`\sigma_u`, and :math:`\sigma_{\epsilon}` are estimated via MLE under
the alternative hypothesis.
The association test is based on the comparison of the two marginal likelihoods
learnt under the null and alternative hypotheses.

We now show how to use :mod:`limix` to perform association tests using
linear mixed models.
The outcome-vector is given by ``y``.
The covariance matrix is defined by the ``kinship`` variable.
We do not provide any covariate.
In that case, the function :func:`limix.qtl.scan` we call will internally add
a covariate of ones to be multiplied by the offset parameter.
Finally, we pass a matrix ``candidates`` of four columns representing four
alternative hypotheses to be tested:

.. doctest::

    >>> from numpy.random import RandomState
    >>> from numpy import dot
    >>> from limix.qtl import scan
    >>>
    >>> random = RandomState(1)
    >>>
    >>> n = 100
    >>>
    >>> candidates = (random.rand(n, 4) < 0.2).astype(float)
    >>> y = random.randn(n)
    >>> X = random.randn(n, 10)
    >>> kinship = dot(X, X.T) / float(10)
    >>>
    >>> model = scan(candidates, y, 'normal', kinship, verbose=False)
    >>> print(model.variant_pvalues)  # doctest: +NPY_FLEX_NUMS
    0    0.128600
    1    0.210929
    2    0.231456
    3    0.037974
    dtype: float64
    >>> print(model.variant_effsizes)  # doctest: +NPY_FLEX_NUMS
    0    0.314872
    1   -0.245231
    2   -0.182707
    3    0.511300
    dtype: float64
    >>> print(model.variant_effsizes_se)  # doctest: +NPY_FLEX_NUMS
    0    0.207201
    1    0.196026
    2    0.152686
    3    0.246394
    dtype: float64
    >>> print(model)  # doctest: +NPY_FLEX_NUMS
    Variants
           effsizes  effsizes_se   pvalues
    count  4.000000     4.000000  4.000000
    mean   0.099559     0.200577  0.152240
    std    0.371683     0.038545  0.088195
    min   -0.245231     0.152686  0.037974
    25%   -0.198338     0.185191  0.105943
    50%    0.066083     0.201614  0.169765
    75%    0.363979     0.216999  0.216061
    max    0.511300     0.246394  0.231456
    <BLANKLINE>
    Covariate effect sizes for the null model
      offset
     0.00139

The above example prints the estimated p-value, effect size, and standard
error of the effect size of each variant.
It also shows a summary of the result by printing the variable ``model``, an
instance of the :class:`limix.qtl.model.QTLModel` class.

A **generalised linear mixed model** (GLMM) in an extension of a LMM that
allows for residual errors distributed according to an exponential-family
distribution.
Let us replace :math:`\mathbf y` in the LMM equation by :math:`\mathbf z`, and
define the outcome-vector as

.. math::

    y_i ~|~ z_i \sim \text{ExpFam}(\mu_i = g(z_i)).

The multivariate Normal distribution :math:`\mathbf z` is
considered a latent (unobserved) variable.
The :math:`\mu_i` variable is the parameter defining the expected value of a
distribution :math:`\text{ExpFam}(\cdot)`.
It is defined via a link function :math:`g(\cdot)`, which converts the interval
of :math:`z_i` (real numbers) to the appropriate interval for :math:`\mu_i`.

The following example applies :func:`limix.qtl.scan` to perform five
likelihood-ratio tests for association with an outcome vector ``y`` having
residual errors that follow a Poisson distribution.
The matrix ``G`` defines both five alternative hypotheses
(the first five columns) and the covariance matrix (the remaining columns).

.. doctest::

    >>> from numpy import dot, exp, sqrt
    >>> from numpy.random import RandomState
    >>> from limix.qtl import scan
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
    >>> model = scan(candidates, y, 'poisson', K, verbose=False)
    >>>
    >>> print(model.variant_pvalues)
    0    0.073623
    1    0.341348
    2    0.595322
    3    0.742450
    4    0.781789
    dtype: float64
    >>> print(model.variant_effsizes)
    0    2.477110
    1   -1.256363
    2   -0.705970
    3   -0.476551
    4    0.376326
    dtype: float64
    >>> print(model.variant_effsizes_se)
    0    1.384672
    1    1.320395
    2    1.329161
    3    1.450197
    4    1.358630
    dtype: float64
    >>> print(model)
    Variants
           effsizes  effsizes_se   pvalues
    count  5.000000     5.000000  5.000000
    mean   0.082910     1.368611  0.506906
    std    1.461785     0.052189  0.297491
    min   -1.256363     1.320395  0.073623
    25%   -0.705970     1.329161  0.341348
    50%   -0.476551     1.358630  0.595322
    75%    0.376326     1.384672  0.742450
    max    2.477110     1.450197  0.781789
    <BLANKLINE>
    Covariate effect sizes for the null model
       offset
     0.007468

Out-of-core
^^^^^^^^^^^

Limix is also able to handle large datasets that do not fit in memory.
This is internally accomplished by mapping chunks of an array to the
corresponding parts of the files that define the dataset.
Those chunks of memory are then instantiated (and then released) on demand
by the underlying algorithm.

Interface
^^^^^^^^^

.. autofunction:: limix.qtl.scan
.. autoclass:: limix.qtl.model.QTLModel
    :members:
.. autofunction:: limix.qtl.iscan
.. autoclass:: limix.qtl.model.IQTLModel
    :members:

"""

# TODO: move this documentation to rst

from .interact import iscan
from .scan import scan

__all__ = ['scan', 'iscan']
