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
    >>> from pandas import option_context
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
    >>> with option_context('precision', 5):
    ...     print(model.variant_pvalues.to_string())
    0    0.34868
    1    1.00000
    2    0.42545
    3    0.05915
    >>> with option_context('precision', 5):
    ...     print(model.variant_effsizes.to_string())
    0    0.11475
    1    0.00489
    2   -0.20052
    3    0.51587
    >>> with option_context('precision', 5):
    ...     print(model.variant_effsizes_se.to_string())
    0    1.22450e-01
    1    3.02424e+13
    2    2.51589e-01
    3    2.73372e-01
    >>> with option_context('precision', 5):
    ...     print(model)
    Variants
           effsizes  effsizes_se  pvalues
    count   4.00000  4.00000e+00  4.00000
    mean    0.10875  7.56059e+12  0.45832
    std     0.30123  1.51212e+13  0.39406
    min    -0.20052  1.22450e-01  0.05915
    25%    -0.04646  2.19304e-01  0.27630
    50%     0.05982  2.62480e-01  0.38706
    75%     0.21503  7.56059e+12  0.56908
    max     0.51587  3.02424e+13  1.00000
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
    >>> with option_context('precision', 5):
    ...     print(model.variant_pvalues.to_string())
    0    0.06563
    1    0.31398
    2    0.54154
    3    0.65972
    4    0.68753
    >>> with option_context('precision', 5):
    ...     print(model.variant_effsizes.to_string())
    0    2.47711
    1   -1.25636
    2   -0.70597
    3   -0.47655
    4    0.37633
    >>> with option_context('precision', 5):
    ...     print(model.variant_effsizes_se.to_string())
    0    1.34557
    1    1.24775
    2    1.15640
    3    1.08233
    4    0.93565
    >>> with option_context('precision', 5):
    ...     print(model)
    Variants
           effsizes  effsizes_se  pvalues
    count   5.00000      5.00000  5.00000
    mean    0.08291      1.15354  0.45368
    std     1.46179      0.15677  0.26218
    min    -1.25636      0.93565  0.06563
    25%    -0.70597      1.08233  0.31398
    50%    -0.47655      1.15640  0.54154
    75%     0.37633      1.24775  0.65972
    max     2.47711      1.34557  0.68753
    <BLANKLINE>
    Covariate effect sizes for the null model
      offset
     0.00747

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
