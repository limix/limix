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
Covariates are defined by the columns of :math:`\mathrm M` and
:math:`\mathrm X` commonly contain all genetic variants of each sample.

The outcome-vector is thus distributed according to

.. math::

    \mathbf y \sim \mathcal N(\mathrm M\boldsymbol\beta_0,
                              \sigma_u^2 \mathrm X \mathrm X^{\intercal}
                              + \sigma_e^2\mathrm I).

The parameters :math:`\boldsymbol\beta_0`, :math:`\sigma_u`, and
:math:`\sigma_{\epsilon}` are estimated via the maximum likelihood estimation
(MLE) approach under the null hypothesis.

The alternative hypothesis for single-variant testing consists in the addition
of a fixed-effect size :math:`\beta_1`:

.. math::

    \mathbf y = \mathrm M\boldsymbol\beta_1 + \mathbf g\beta_1
        + \mathrm X\mathbf u + \boldsymbol\epsilon.

The parameter :math:`\beta_1` multiplies a given vector :math:`\mathbf g`,
typically representing a genetic locus of interest.

The following example applies :func:`limix.qtl.scan` to perform four
likelihood-ratio tests for association with an outcome vector ``y``:

.. doctest::

    >>> from numpy.random import RandomState
    >>> from numpy import dot
    >>> from limix.qtl import scan
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
    >>> print(model.variant_pvalues)
    [ 0.3487  1.      0.4254  0.0592]
    >>> print(model.variant_effsizes)
    [ 0.1148  0.0049 -0.2005  0.5159]
    >>> print(model.variant_effsizes_se)
    [  1.2245e-01   3.0242e+13   2.5159e-01   2.7337e-01]
    >>> print(model)
    Variants
           effsizes   effsizes_se   pvalues
    count  4.000000  4.000000e+00  4.000000
    mean   0.108750  7.560590e+12  0.458319
    std    0.301228  1.512118e+13  0.394062
    min   -0.200519  1.224501e-01  0.059153
    25%   -0.046459  2.193043e-01  0.276297
    50%    0.059825  2.624805e-01  0.387062
    75%    0.215033  7.560590e+12  0.569085
    max    0.515868  3.024236e+13  1.000000
    <BLANKLINE>
    Covariate effect sizes for the null model
        offset
    0.00139039

The above example prints the estimated p-values, effect sizes, and standard
errors of the effect sizes of each variant.
It also shows a summary of the result by printing the variable ``model``, an
instance of :class:`limix.qtl.model.QTLModel`.

A **generalised linear mixed model** (GLMM) in an extension of a LMM that allows
for residual errors distributed according to an exponential-family
distribution.
Let us replace :math:`\mathbf y` in the LMM equation by :math:`\mathbf z`, and
define the outcome-vector as

.. math::

    y_i ~|~ z_i \sim \text{ExpFam}(\mu_i = g(z_i)).

In other words, the multivariate Normal distribution :math:`\mathbf z` is
considered a latent (unobserved) variable.
It defines the expected value of the outcome via a (link) function
:math:`g(\cdot)`.

The following example applies :func:`limix.qtl.scan` to perform five
likelihood-ratio tests for association with an outcome vector ``y`` having
residual errors that follow a Poisson distribution:

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
    [ 0.0694  0.3336  0.5899  0.7387  0.7796]
    >>> print(model.variant_effsizes)
    [ 2.4732 -1.2588 -0.7068 -0.4772  0.3752]
    >>> print(model.variant_effsizes_se)
    [ 1.362   1.3018  1.3112  1.4309  1.3405]
    >>> print(model)
    Variants
           effsizes  effsizes_se   pvalues
    count  5.000000     5.000000  5.000000
    mean   0.081131     1.349272  0.502224
    std    1.460864     0.051503  0.298472
    min   -1.258785     1.301777  0.069380
    25%   -0.706753     1.311199  0.333557
    50%   -0.477232     1.340503  0.589878
    75%    0.375201     1.361954  0.738747
    max    2.473225     1.430925  0.779557
    <BLANKLINE>
    Covariate effect sizes for the null model
        offset
    0.00746817

Interface
^^^^^^^^^

.. autofunction:: limix.qtl.scan
.. autoclass:: limix.qtl.model.QTLModel
    :members:

"""

from .scan import scan

__all__ = ['scan']
