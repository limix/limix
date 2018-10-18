************************
Quantitative trait locus
************************


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
    >>> print(model.variant_pvalues.to_dataframe()) # doctest: +FLOAT_CMP
                     pv
    candidate
    0           0.74981
    1           0.00537
    2           0.07036
    3           0.97155
    >>> print(model.variant_effsizes.to_dataframe()) # doctest: +FLOAT_CMP
                effsizes
    candidate
    0          -0.09660
    1          -1.02874
    2          -0.46314
    3          -0.01155
    >>> print(model.variant_effsizes_se.to_dataframe()) # doctest: +FLOAT_CMP
               effsizes std
    candidate
    0               0.30293
    1               0.36956
    2               0.25594
    3               0.32377
    >>> print(model) # doctest: +FLOAT_CMP
    Variants
    --------
           effsizes  effsizes_se  pvalues
    count   4.00000      4.00000  4.00000
    mean   -0.40001      0.31305  0.44927
    std     0.46269      0.04716  0.48433
    min    -1.02874      0.25594  0.00537
    25%    -0.60454      0.29118  0.05411
    50%    -0.27987      0.31335  0.41008
    75%    -0.07534      0.33522  0.80524
    max    -0.01155      0.36956  0.97155
    <BLANKLINE>
    Covariate effect sizes for H0
    -----------------------------
     0
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
    >>> print(model.variant_pvalues.to_dataframe()) # doctest: +FLOAT_CMP
                    pv
    candidate
    0          0.19819
    1          0.44134
    2          0.47341
    3          0.21548
    4          0.70666
    >>> print(model.variant_effsizes.to_dataframe()) # doctest: +FLOAT_CMP
               effsizes
    candidate
    0           1.69168
    1          -1.00863
    2          -1.24902
    3           2.04198
    4          -0.50974
    >>> print(model.variant_effsizes_se.to_dataframe()) # doctest: +FLOAT_CMP
               effsizes std
    candidate
    0               1.31470
    1               1.31003
    2               1.74216
    3               1.64859
    4               1.3544
    >>> print(model) # doctest: +FLOAT_CMP
    Variants
    --------
           effsizes  effsizes_se  pvalues
    count   5.00000      5.00000  5.00000
    mean    0.19325      1.47399  0.40702
    std     1.55579      0.20552  0.20956
    min    -1.24902      1.31003  0.19819
    25%    -1.00863      1.31470  0.21548
    50%    -0.50974      1.35444  0.44134
    75%     1.69168      1.64859  0.47341
    max     2.04198      1.74216  0.70666
    <BLANKLINE>
    Covariate effect sizes for H0
    -----------------------------
     0
       -0.00042

Interface
^^^^^^^^^

.. autofunction:: limix.qtl.scan
        :noindex:
.. autoclass:: limix.qtl.QTLModel
    :members:
    :noindex:

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
