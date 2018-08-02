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
The comparison of the two marginal likelihoods learnt under the null and alternative
hypotheses allows us to perform a likelihood ratio test [LRT]_.

We now show how to use limix to perform association tests using
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
    >>> n = 25
    >>>
    >>> candidates = (random.rand(n, 4) < 0.2).astype(float)
    >>> y = random.randn(n)
    >>> X = random.randn(n, 50)
    >>> kinship = dot(X, X.T) / float(25)
    >>>
    >>> model = scan(candidates, y, 'normal', kinship, verbose=False)
    >>> print(model.variant_pvalues) # doctest: +FLOAT_CMP
    candidate0    0.749809
    candidate1    0.005375
    candidate2    0.070358
    candidate3    0.971545
    dtype: float64
    >>> print(model.variant_effsizes) # doctest: +FLOAT_CMP
    candidate0   -0.096603
    candidate1   -1.028743
    candidate2   -0.463141
    candidate3   -0.011549
    dtype: float64
    >>> print(model.variant_effsizes_se) # doctest: +FLOAT_CMP
    candidate0    0.302935
    candidate1    0.369561
    candidate2    0.255936
    candidate3    0.323771
    dtype: float64
    >>> print(model) # doctest: +FLOAT_CMP
    Variants
           effsizes  effsizes_se   pvalues
    count  4.000000     4.000000  4.000000
    mean  -0.400009     0.313051  0.449272
    std    0.462691     0.047162  0.484325
    min   -1.028743     0.255936  0.005375
    25%   -0.604541     0.291185  0.054112
    50%   -0.279872     0.313353  0.410084
    75%   -0.075340     0.335218  0.805243
    max   -0.011549     0.369561  0.971545
    <BLANKLINE>
    Covariate effect sizes for the null model
     covariate0
        0.04828

The above example prints the estimated p-value, effect size, and standard
error of the effect size of each variant.
It also shows a summary of the result by printing the variable ``model``, an
instance of the :class:`limix.qtl.QTLModel` class.

A **generalised linear mixed model** (GLMM) [McC89]_ [McC11]_ in an extension of a LMM
that allows for residual errors distributed according to an exponential-family
distribution [ExFam]_.
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
likelihood ratio tests for association with an outcome vector ``y`` having
residual errors that follow a Poisson distribution.
The matrix ``G`` defines both the five alternative hypotheses
(the first five columns) and the covariance matrix (the remaining columns).

.. doctest::

    >>> from numpy import dot, exp, sqrt
    >>> from numpy.random import RandomState
    >>> from limix.qtl import scan
    >>>
    >>> random = RandomState(0)
    >>>
    >>> G = random.randn(25, 50) / sqrt(50)
    >>> beta = 0.01 * random.randn(50)
    >>>
    >>> z = dot(G, beta) + 0.1 * random.randn(25)
    >>> z += dot(G[:, 0], 1) # causal SNP
    >>>
    >>> y = random.poisson(exp(z))
    >>>
    >>> candidates = G[:, :5]
    >>> K = dot(G[:, 5:], G[:, 5:].T)
    >>> model = scan(candidates, y, 'poisson', K, verbose=False)
    >>>
    >>> print(model.variant_pvalues) # doctest: +FLOAT_CMP
    candidate0    0.198186
    candidate1    0.441343
    candidate2    0.473412
    candidate3    0.215485
    candidate4    0.706657
    dtype: float64
    >>> print(model.variant_effsizes) # doctest: +FLOAT_CMP
    candidate0    1.691677
    candidate1   -1.008630
    candidate2   -1.249025
    candidate3    2.041978
    candidate4   -0.509744
    dtype: float64
    >>> print(model.variant_effsizes_se) # doctest: +FLOAT_CMP
    candidate0    1.314702
    candidate1    1.310032
    candidate2    1.742165
    candidate3    1.648586
    candidate4    1.354442
    dtype: float64
    >>> print(model) # doctest: +FLOAT_CMP
    Variants
           effsizes  effsizes_se   pvalues
    count  5.000000     5.000000  5.000000
    mean   0.193251     1.473985  0.407016
    std    1.555792     0.205517  0.209562
    min   -1.249025     1.310032  0.198186
    25%   -1.008630     1.314702  0.215485
    50%   -0.509744     1.354442  0.441343
    75%    1.691677     1.648586  0.473412
    max    2.041978     1.742165  0.706657
    <BLANKLINE>
    Covariate effect sizes for the null model
     covariate0
       -0.00042

Interface
^^^^^^^^^

.. autofunction:: limix.qtl.scan
.. autoclass:: limix.qtl.QTLModel
    :members:

.. rubric:: References

.. [McC89] McCullagh, Peter, and John A. Nelder. Generalized linear models. Vol. 37.
           CRC press, 1989.
.. [McC11] McCulloch, Charles E., and Shayle R. Searle. Generalized, linear, and mixed
           models. John Wiley & Sons, 2004.
.. [ExFam] Wikipedia contributors. (2018, June 29). Exponential family. In Wikipedia,
           The Free Encyclopedia. Retrieved 13:47, July 26, 2018, from
           https://en.wikipedia.org/w/index.php?title=Exponential_family&oldid=848114709
.. [LRT]   Wikipedia contributors. (2018, June 6). Likelihood-ratio test. In Wikipedia,
           The Free Encyclopedia. Retrieved 13:50, July 26, 2018, from
           https://en.wikipedia.org/w/index.php?title=Likelihood-ratio_test&oldid=844734768
