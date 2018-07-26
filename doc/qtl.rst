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
    >>> n = 100
    >>>
    >>> candidates = (random.rand(n, 4) < 0.2).astype(float)
    >>> y = random.randn(n)
    >>> X = random.randn(n, 10)
    >>> kinship = dot(X, X.T) / float(10)
    >>>
    >>> model = scan(candidates, y, 'normal', kinship, verbose=False)
    >>> print(model.variant_pvalues) # doctest: +FLOAT_CMP
    candidate0    0.227204
    candidate1    0.398132
    candidate2    0.446016
    candidate3    0.063060
    dtype: float64
    >>> print(model.variant_effsizes) # doctest: +FLOAT_CMP
    candidate0    0.314881
    candidate1   -0.245241
    candidate2   -0.182706
    candidate3    0.511296
    dtype: float64
    >>> print(model.variant_effsizes_se) # doctest: +FLOAT_CMP
    candidate0    0.260751
    candidate1    0.290239
    candidate2    0.239749
    candidate3    0.275073
    dtype: float64
    >>> print(model) # doctest: +FLOAT_CMP
    Variants
           effsizes  effsizes_se   pvalues
    count  4.000000     4.000000  4.000000
    mean   0.099557     0.266453  0.283603
    std    0.371686     0.021492  0.174466
    min   -0.245241     0.239749  0.063060
    25%   -0.198340     0.255500  0.186168
    50%    0.066087     0.267912  0.312668
    75%    0.363985     0.278864  0.410103
    max    0.511296     0.290239  0.446016
    <BLANKLINE>
    Covariate effect sizes for the null model
     covariate0
      -0.097162

The above example prints the estimated p-value, effect size, and standard
error of the effect size of each variant.
It also shows a summary of the result by printing the variable ``model``, an
instance of the :class:`limix.qtl.model.QTLModel` class.

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
    >>> print(model.variant_pvalues) # doctest: +FLOAT_CMP
    candidate0    0.073623
    candidate1    0.341348
    candidate2    0.595322
    candidate3    0.742450
    candidate4    0.781789
    dtype: float64
    >>> print(model.variant_effsizes) # doctest: +FLOAT_CMP
    candidate0    2.477110
    candidate1   -1.256363
    candidate2   -0.705970
    candidate3   -0.476551
    candidate4    0.376326
    dtype: float64
    >>> print(model.variant_effsizes_se) # doctest: +FLOAT_CMP
    candidate0    1.384672
    candidate1    1.320395
    candidate2    1.329161
    candidate3    1.450197
    candidate4    1.358630
    dtype: float64
    >>> print(model) # doctest: +FLOAT_CMP
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
     covariate0
       0.007468

Interface
^^^^^^^^^

.. autofunction:: limix.qtl.scan
.. autoclass:: limix.qtl.model.QTLModel
    :members:
.. autofunction:: limix.qtl.iscan
.. autoclass:: limix.qtl.model.IQTLModel
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
