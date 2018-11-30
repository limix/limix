************************
Quantitative trait locus
************************


Introduction
^^^^^^^^^^^^

Every genetic model considered here is an instance of **generalised linear mixed model**
(GLMM).
It consists in four main components [St16]_:

    - A linear predictor, :math:`\mathbf z = \mathbf M\boldsymbol\beta + \mathbf
      X\mathbf u`.
    - The distribution of the random effects, :math:`\mathbf u \sim \mathcal N(\mathbf
      0, \boldsymbol\Sigma)`.
    - The distribution of the outcome conditioned on the random effects (also known as
      the residual distribution), :math:`y_i | \mathbf u`.
    - The link function, :math:`g(\mu_i) = z_i`.

The term :math:`\mu_i` represents the mean of :math:`y_i` conditioned on :math:`\mathbf
u`:

.. math::

    \mu_i = \mathbb E[y_i | \mathbf u].

The role of the link function is to scale the domain of :math:`z_i`, which ranges from
:math:`-\infty` to :math:`+infty`, to the redisual distrubution parameter :math:`\mu_i`.
For example, the mean of a Bernoulli distribution is bounded within :math:`[0, 1]`, and
therefore requires a link function to translate values of :math:`z_i` into values of
:math:`\mu_i`.

The distribution of the outcome, conditioned on the random effects, has to be one from
the exponential family [Ef18]_ having mean :math:`\mu_i`:

.. math::

    y_i | \mathbf u \sim \text{ExpFam}(\mu_i).

A notable instance of the above model is the **linear mixed model** (LMM).
It consists of the identity link function, :math:`g(\mu_i) = \mu_i`, and of normally
distributed residuals, :math:`y_i | \mathbf u \sim \mathcal N(\mu_i, \sigma_i^2)`
[Mc11]_.
It is more commonly described by the equation

.. math::
    :label: eq_lmm

    \mathbf y = \mathbf M\boldsymbol\beta + \mathbf X\mathbf u + \boldsymbol\epsilon,

for which :math:`\epsilon_i\sim\mathcal N(0, \sigma_i^2)`.  The random variables
:math:`\mathbf u` and :math:`\boldsymbol\epsilon` are independent from each other as
well as :math:`\epsilon_i` and :math:`\epsilon_j` for :math:`i\neq j`.  Defining
:math:`\mathbf v = \mathbf X\mathbf u` leads to

.. math::

    \mathbf v \sim \mathcal(\mathbf 0, \mathbf X\boldsymbol\Sigma\mathbf X^{\intercal}).

There is another even simpler instance of GLMM that is also used in genetic analysis:
a **linear model** (LM) is merely a LMM without the random effects:

.. math::

    \mathbf y = \mathbf M\boldsymbol\beta + \boldsymbol\epsilon.

The above models are used to establish a statiscal tests to find significant association
between genetic loci and phenotype.
For that, their parameters have to be estimated.

As an example, let us define two parameters that will describe the overall variances of
the random effects
and of the residual effects:

.. math::

    \boldsymbol\Sigma = v_0\mathbf I_0 ~~\text{and}~~
    \sigma_i^2 = v_1.

If we assume a LMM, this example of model can be described by Eq. :eq:`eq_lmm` for which

.. math::

    \mathbf v\sim\mathcal N(\mathbf 0, v_0\mathbf X\mathbf X^{\intercal}) ~~\text{and}~~
    \boldsymbol\epsilon\sim\mathcal N(\mathbf 0, v_1\mathbf I_1).

Therefore we have a model with three parameters: an array of effect sizes
:math:`\boldsymbol\beta`
and variances :math:`v_0` and :math:`v_1`.

Statistical test
^^^^^^^^^^^^^^^^

We use the **likelihood ratio test** (LRT) approach [LR18]_ to assess the significance
of the association
between genetic variants and the phenotype.
It is based on the ratio between the marginal likelihood of the null and alternative
models:

.. math::

    \mathcal H_0: \boldsymbol\theta_0\\
    \mathcal H_1: \boldsymbol\theta_1

where :math:`\boldsymbol\theta_0` is fit under the null model and
:math:`\boldsymbol\theta_1` is fit under the alternative model.
The parameter inference is done via the maximum likelihood estimation (MLE) approach
[ML18]_:

.. math::

    \boldsymbol{\hat{\theta}} = \underset{\boldsymbol\theta}{\mathrm{argmax}}~~
        p(\mathbf y | \mathbf M, \mathbf X; \boldsymbol\theta).

The likelihood ratio is then equal to

.. math::

    \frac{p(\mathbf y| \mathbf M, \mathbf X; \boldsymbol{\hat{\theta_0}})}
        {p(\mathbf y| \mathbf M, \mathbf X; \boldsymbol{\hat{\theta_1}})}.

which will define the p-value of that comparison.

Single-trait association
^^^^^^^^^^^^^^^^^^^^^^^^

We first consider that the observed phenotype is described by additive effects from
covariates and genetic components, and any deviation from that is captured by the
assumed residual distribution and/or an over-dispersion component.  Let :math:`\mathbf
M` be a matrix of covariates and let :math:`\mathbf G` be a matrix of genetic variants
that we suspect might have some effect on the phenotype.  Therefore, we have the linear
model

.. math::

    \mathbf y = \underbrace{\mathbf M\boldsymbol\alpha}_{\text{covariates}}+
    \underbrace{\mathbf G\boldsymbol\beta}_{\text{genetics}}+
    \underbrace{\boldsymbol\epsilon}_{\text{noise}},\\
    \text{where}~~\boldsymbol\epsilon\sim\mathcal N(\mathbf 0, v_1\mathbf I),~~~~~~

and we wish to compare the following hypotheses:

.. math::

    \mathcal H_0: \boldsymbol\beta = 0\\
    \mathcal H_1: \boldsymbol\beta \neq 0

Note that the parameters of the above model are the covariate effect sizes,
:math:`\boldsymbol\alpha`, the effect sizes of a set of genetic variants,
:math:`\boldsymbol\beta`, and the variance :math:`v_1` of the noise variable.  Under the
null hypothesis, we set :math:`\boldsymbol\beta=\mathbf 0` and fit the rest of the
parameters.  Under the alternative hypothesis, we learn all the parameters.  At the end,
we compare the marginal likelihoods via the likelihood ratio test.

Let us first generate a random dataset having a phenotype, covariates, and a set of
genetic candidates.

.. doctest::

    >>> from numpy.random import RandomState
    >>> from numpy import dot, ones, stack
    >>> from pandas import DataFrame
    >>>
    >>> random = RandomState(1)
    >>>
    >>> # 25 samples
    >>> n = 25
    >>>
    >>> # genetic variants
    >>> snps = (random.rand(n, 4) < 0.2).astype(float)
    >>>
    >>> #phenotype
    >>> y = random.randn(n)
    >>>
    >>> # offset
    >>> offset = ones(n)
    >>> # age
    >>> age = random.randint(16, 75, n)
    >>> M = DataFrame(stack([offset, age], axis=1), columns=["offset", "age"])
    >>> print(M.head())
        offset      age
    0  1.00000 49.00000
    1  1.00000 18.00000
    2  1.00000 36.00000
    3  1.00000 35.00000
    4  1.00000 64.00000

We now apply the function :func:`limix.qtl.st_scan` to our dataset

.. doctest::

    >>> from limix.qtl import st_scan
    >>>
    >>> r = st_scan(snps, y, 'normal', M=M, verbose=False)
    >>> print(r)
    Variants
    --------
            effsizes  effsizes_se  pvalues
    count         4            4        4
    mean   -0.28530      0.28270  0.51540
    std     0.37510      0.04897  0.49548
    min    -0.79864      0.23277  0.02219
    25%    -0.44725      0.25768  0.12263
    50%    -0.17812      0.27441  0.53896
    75%    -0.01616      0.29944  0.93172
    max     0.01366      0.34920  0.96147
    <BLANKLINE>
    Covariate effect sizes for H0
    -----------------------------
         age   offset
     0.02016 -0.81859

The variable ``r`` is instance of the class :class:`limix.qtl.QTLResult` and stores all
the results of the analysis.  Printing it as we did above it will show a summary of the
results.

Suppose we also have access to the whole genotype of our samples, :math:`\mathbf X`, and
we want to use them to account for population structure and cryptic relatedness in our
data (CITE).  Since the number of genetic variants in :math:`\mathbf X` is commonly
larger than the number of samples, and because we are not acctually interested in their
effect sizes, we will include it in our model as a random component.  We now have a
**linear mixed model**:

.. math::

    \mathbf y = \underbrace{\mathbf M\boldsymbol\alpha}_{\text{covariates}}+
    \underbrace{\mathbf G\boldsymbol\beta}_{\text{genetics}}+
    \underbrace{\mathbf X\mathbf u}_{\text{pop. struct.}}+
    \underbrace{\boldsymbol\epsilon}_{\text{noise}},\\
    \text{where}~~
        \mathbf u\sim\mathcal N(\mathbf 0, v_0\mathbf I_0) ~~\text{and}
    ~~\boldsymbol\epsilon\sim\mathcal N(\mathbf 0, v_1\mathbf I_1).

It is important to note that :math:`\mathbf v=\mathbf X\mathbf u` can be equivalenty
described by a multivariate Normal distribution with a covariance proportional to
:math:`\mathbf K = \mathbf X\mathbf X^{\intercal}`:

.. math::

    \mathbf v \sim \mathcal N(\mathbf 0, v_0\mathbf K).

We perform the analysis again now using also the covariance :math:`\mathbf K` by calling
the function :func:`limix.qtl.st_scan`.

.. doctest::

    >>> from limix.stats import linear_kinship
    >>>
    >>> # Whole genotype of each sample.
    >>> X = random.randn(n, 50)
    >>> # Estimate a kinship relationship betweem samples.
    >>> K = linear_kinship(X, verbose=False)
    >>>
    >>> model = st_scan(X, y, 'normal', K, M=M, verbose=False)
    >>> print(model.variant_pvalues.to_dataframe().head()) # doctest: +FLOAT_CMP
                    pv
    candidate
    0          0.72152
    1          0.61718
    2          0.44223
    3          0.94542
    4          0.40416
    >>> print(model.variant_effsizes.to_dataframe().head()) # doctest: +FLOAT_CMP
               effsizes
    candidate
    0           0.04675
    1          -0.05855
    2          -0.09668
    3           0.00746
    4           0.12734
    >>> print(model.variant_effsizes_se.to_dataframe().head()) # doctest: +FLOAT_CMP
               effsizes std
    candidate
    0               0.13116
    1               0.11713
    2               0.12582
    3               0.10899
    4               0.15264
    >>> print(model) # doctest: +FLOAT_CMP
    Variants
    --------
           effsizes  effsizes_se  pvalues
    count        50           50       50
    mean    0.00412      0.12327  0.51505
    std     0.12437      0.01564  0.29979
    min    -0.25728      0.09544  0.01043
    25%    -0.08632      0.11008  0.31724
    50%     0.00039      0.12209  0.50611
    75%     0.07436      0.13173  0.76775
    max     0.31626      0.16753  0.97555
    <BLANKLINE>
    Covariate effect sizes for H0
    -----------------------------
         age   offset
     0.02016 -0.81859


Generalised phenotype
~~~~~~~~~~~~~~~~~~~~~

If the residuals of the phenotype does not follow a Normal distribution, then we might
consider perform analysis using a **generalised linear mixed model**.  Let us consider
Poisson distributed residuals:

.. math::

    y_i | \mathbf z \sim \text{Bernoulli}(g(\mu_i)=z_i).

In the latter case, the :math:`\boldsymbol\epsilon` can be used to describe the
dispersion between samples not fully captured by the residual distribution.

The following example applies :func:`limix.qtl.st_scan` to perform five likelihood ratio
tests for association with an outcome vector ``y`` having residual errors that follow a
Poisson distribution.  The matrix ``G`` defines both the five alternative hypotheses
(the first five columns) and the covariance matrix (the remaining columns).

.. doctest::

    >>> from numpy import exp, sqrt
    >>> from numpy.random import RandomState
    >>> from limix.qtl import st_scan
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
    >>> K = linear_kinship(G[:, 5:], verbose=False)
    >>> model = st_scan(candidates, y, 'poisson', K, verbose=False)
    >>>
    >>> print(model.variant_pvalues.to_dataframe()) # doctest: +FLOAT_CMP
                    pv
    candidate
    0          0.21645
    1          0.44194
    2          0.48394
    3          0.23413
    4          0.69534
    >>> print(model.variant_effsizes.to_dataframe()) # doctest: +FLOAT_CMP
               effsizes
    candidate
    0           1.62727
    1          -1.02366
    2          -1.23573
    3           1.97540
    4          -0.53729
    >>> print(model.variant_effsizes_se.to_dataframe()) # doctest: +FLOAT_CMP
               effsizes std
    candidate
    0               1.31655
    1               1.33129
    2               1.76537
    3               1.66030
    4               1.37198
    >>> print(model) # doctest: +FLOAT_CMP
    Variants
    --------
           effsizes  effsizes_se  pvalues
    count         5            5        5
    mean    0.16120      1.48910  0.41436
    std     1.52348      0.20859  0.19761
    min    -1.23573      1.31655  0.21645
    25%    -1.02366      1.33129  0.23413
    50%    -0.53729      1.37198  0.44194
    75%     1.62727      1.66030  0.48394
    max     1.97540      1.76537  0.69534
    <BLANKLINE>
    Covariate effect sizes for H0
    -----------------------------
      offset
    -0.01412

Single-trait with interaction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following linear mixed model is considered:

.. math::

    \mathbf{y} =
    \underbrace{\mathbf M\boldsymbol\beta}_
            {\substack{\text{fixed effects}\\ \text{without interaction}}}+
    \underbrace{(\mathbf G\odot\mathbf E_0)\boldsymbol\beta_0}_{\mathrm G\times\mathrm E_0} +
    \underbrace{\mathbf G\odot\mathbf E_1\boldsymbol\beta_1}_{\mathrm G\times\mathrm E_1} +
    \underbrace{\mathbf X\mathbf{u}}_{\text{random effects}}+
    \underbrace{\boldsymbol{\epsilon}}_{\text{residual}}.

The **GxE** terms are also fixed effects but encoding the interations between genetic
variants and environmental covariates defined by the user.

.. doctest::

    >>> from numpy import concatenate, newaxis
    >>> from limix.qtl import st_iscan
    >>> # generate interacting variables (environment)
    >>> random = RandomState(1)
    >>> E = random.randn(y.shape[0], 1)
    >>>
    >>> # add additive environment as covariate
    >>> ME = concatenate([M, E], axis=1)
    >>>
    >>> snps = random.randn(n, 100)
    >>>
    >>> # interaction test
    >>> res = st_iscan(snps, y[:, newaxis], M=ME, E1=E, verbose=False)
    >>> print(res.head())  # doctest: +FLOAT_CMP
           pv1      pv0       pv    beta0  beta0_ste     lrt1     lrt0      lrt
    0  0.14584  0.06186  0.54644  0.36731    0.19671  3.85044  3.48671  0.36373
    1  0.81134  0.52514  0.90466  0.13994    0.22022  0.41813  0.40378  0.01435
    2  0.74902  0.45885  0.86414  0.17079    0.23056  0.57798  0.54871  0.02928
    3  0.77165  0.79650  0.50141 -0.05915    0.22937  0.51845  0.06650  0.45195
    4  0.81176  0.64857  0.64725  0.09675    0.21229  0.41709  0.20770  0.20939


The process method returns three sets of P values: (i) ``pv0`` are association test P
values (:math:`\boldsymbol{\alpha}\neq{0}` when :math:`\boldsymbol{\beta}={0}`), (ii)
``pv1`` are association + interaction P values (:math:`\left[\boldsymbol{\beta},
\boldsymbol{\alpha}\right]\neq{0}`) and (iii) ``pv`` are interaction P values
(:math:`\boldsymbol{\alpha}\neq{0}`).  The effect sizes of the association test are also
returned.

If ``E0`` is not specified, a column-vector of ones is considered.  In this case the
:math:`\mathbf G\odot\mathbf E_0` term reduces to an additive genetic effect, and thus
the test corresponds to a standard gxe test.

If iter0 is provided,

.. doctest::

    >>> # generate interacting variables to condition on
    >>> E0 = random.randn(y.shape[0], 1)
    >>>
    >>> # generate interacting variables to test
    >>> E1 = random.randn(y.shape[0], 1)
    >>>
    >>> # add additive environment as covariate
    >>> ME = concatenate([M, E0, E1], axis=1)
    >>>
    >>> # interaction test
    >>> r = st_iscan(snps, y[:, newaxis], M=ME, E1=E1, E0=E0, verbose=False)
    >>> print(r.head())  # doctest: +FLOAT_CMP
           pv1      pv0       pv     lrt1     lrt0      lrt
    0  0.36534  0.22031  0.47451  2.01383  1.50237  0.51146
    1  0.28558  0.15232  0.49876  2.50648  2.04891  0.45757
    2  0.18256  0.07701  0.60042  3.40136  3.12700  0.27436
    3  0.61833  0.57460  0.42139  0.96148  0.31504  0.64644
    4  0.71350  0.47786  0.67886  0.67515  0.50374  0.17142



StructLMM
^^^^^^^^^

StructLMM can be use to test for interaction with multiple environments or to test for
association of genetic variants while accounting for GxE interactions.
The StructLMM model is

.. math::
    \mathbf{y}=
    \underbrace{\mathbf{M}\boldsymbol\beta}_{\text{covariates}}+
    \underbrace{\mathbf{x}\odot\boldsymbol\gamma}_{\text{genetics}}+
    \underbrace{\mathbf E\mathbf u}_{\text{random effects}}+
    \underbrace{\boldsymbol\epsilon}_{\text{noise}},

where

.. math::
    \boldsymbol\gamma\sim\mathcal N(\mathbf 0,
    \sigma^2_g(\underbrace{(1-\rho)\mathbf 1}_{\text{persistent}}
        + \underbrace{\rho\mathbf E\mathbf E^{\intercal}}_{\text{GxE}}),\\
    \mathbf u\sim\mathcal N(\mathbf 0, v_0\mathbf I),
    ~~\text{and}~~
    \boldsymbol\epsilon\sim\mathcal N(\mathbf 0, v_1\mathbf I).

.. doctest::

    >>> from limix.qtl import st_sscan
    >>>
    >>> E = random.randn(y.shape[0], 10)
    >>>
    >>> r = st_sscan(snps[:, :5], y[:, newaxis], E, tests=['inter', 'assoc'],
    ...              verbose=False)
    >>> print(r.head())  # doctest: +FLOAT_CMP
           pvi      pva
    0  0.05753  0.05925
    1  0.14415  0.17995
    2  0.22886  0.34244
    3  0.69552  0.82283
    4  0.59211  0.82515


The process method returns two sets of P values:
(i) ``pvi`` are the interaction P values,
(ii) ``pva`` are the association P values.


Multi-trait association
^^^^^^^^^^^^^^^^^^^^^^^

A multi-trait linear mixed model is still a LMM.
Therefore, its equation

.. math::

    \text{vec}(\mathbf{Y}) =
    \underbrace{(\mathbf A_c \otimes \mathbf M) \text{vec}(\mathbf B_c)}_{\text{covariates}}+
    \underbrace{(\mathbf A_g \otimes \mathbf G) \text{vec}(\mathbf B_g)}_{\text{genetics}}+
    \underbrace{\text{vec}(\mathbf U)}_{\text{random effect}}+
    \underbrace{\text{vec}(\boldsymbol\Psi)}_{\text{noise}}

is equivalent to Eq. :eq:`eq_lmm` but structured in a different way.
The columns of :math:`\mathbf Y` correspond to the different traits being
considered.
The columns of :math:`\mathbf Y` are stacked over each other and is denoted by
:math:`\text{vec}(\mathbf Y)`.
This is a linear transformation called vectorization [Ve19]_, and helps us
describe the model in a more concise manner.

The matrices :math:`\mathbf A_c` and :math:`\mathbf A_g` are design matrices for
the covariates and genetic variants, respectively.
The random effect component is defined by

.. math::

    \text{vec}(\mathbf U)\sim\mathcal N(\mathbf 0, \mathbf C_0\otimes\mathbf K)

and the residuals by

.. math::

    \text{vec}(\boldsymbol\Psi)\sim\mathcal N(\mathbf 0, \mathbf C_1\otimes\mathbf I_n).

As before, :math:`\mathbf M` is the covariates matrix and :math:`\mathbf G` is the
matrix of genetic variants.
The matrices :math:`\mathbf C_0` and :math:`\mathbf C_1` are two matrix-parameters and,
us such, are fitted during the likelihood maximisation.

Any-effect association test
~~~~~~~~~~~~~~~~~~~~~~~~~~~

An any-effect association test corresponds to testing
:math:`\boldsymbol\beta\neq\mathbf 0` with :math:`\mathbf A_g = \mathbf I`.

.. doctest::

    >>> from limix.qtl import mt_scan
    >>> from numpy import eye
    >>>
    >>> p = 4
    >>> Y = random.randn(y.shape[0], p)
    >>>
    >>> Asnps = eye(p)
    >>> r = mt_scan(G, Y, M=M, Asnps=Asnps, K=K, verbose=False)
    >>> print(r.head())  # doctest: +FLOAT_CMP
            pv      lrt
    0  0.79718  1.66438
    1  0.16840  6.44334
    2  0.30726  4.81090
    3  0.46639  3.57615
    4  0.40613  3.99906


Common and interaction tests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The module allows for testing specific trait design matrices for the variant effects.
This is achieved by specifying the two trait design to compare, namely ``Asnps`` and
``Asnps0``.

In the example below we instantiate this principle to test for departures from
a same effect model (same effect size for all analyzed traits).

In this example, the choices of ``Asnps`` and ``Asnps0``
are ``sp.eye(P)`` and ``sp.ones([P, 1])``, respectively.

.. doctest::

    >>> Asnps0 = eye(p)
    >>> r = mt_scan(G, Y, K=K, Ac=None, Asnps=Asnps, Asnps0=Asnps0, verbose=False)
    >>> print(r.head())  # doctest: +FLOAT_CMP
           pv1      pv0       pv     lrt1     lrt0      lrt
    0  0.79947  0.79947      nan  1.65169  1.65169  0.00000
    1  0.15318  0.15318      nan  6.69035  6.69035  0.00000
    2  0.27312  0.27312      nan  5.14113  5.14113  0.00000
    3  0.41205  0.41205      nan  3.95560  3.95560  0.00000
    4  0.39952  0.39952      nan  4.04825  4.04825  0.00000

The process method returns three sets of P values:
(i) ``pv0`` are P values for the association test with snp trait design `Asnps0`,
(ii) ``pv1`` are P values for the association test with snp trait design `Asnps1`,
(iii) ``pv`` are P values for the test `Asnps1` vs `Asnps0`.

In the specific example, these are the P values for
a same-effect association test,
an any-effect association test,
and an any-vs-same effect test.

.. rubric:: References

.. [LR18]  Wikipedia contributors. (2018, October 21). Likelihood-ratio test.
           In Wikipedia, The Free Encyclopedia. Retrieved 16:13, November 27, 2018, from
           https://en.wikipedia.org/w/index.php?title=Likelihood-ratio_test&oldid=865020904
.. [ML18]  Wikipedia contributors. (2018, November 8). Maximum likelihood estimation.
           In Wikipedia, The Free Encyclopedia. Retrieved 16:08, November 27, 2018, from
           https://en.wikipedia.org/w/index.php?title=Maximum_likelihood_estimation&oldid=867823508
.. [St16]  Stroup, W. W. (2016). Generalized linear mixed models: modern concepts, methods
           and applications. CRC press.
.. [Ef18]  Wikipedia contributors. (2018, October 18). Exponential family. In Wikipedia,
           The Free Encyclopedia. Retrieved 18:45, November 25, 2018, from
           https://en.wikipedia.org/w/index.php?title=Exponential_family&oldid=864576150
.. [Mc11]  McCulloch, Charles E., and Shayle R. Searle. Generalized, linear, and mixed
           models. John Wiley & Sons, 2004.
.. [Ve19]  Wikipedia contributors. (2018, September 11). Vectorization (mathematics).
           In Wikipedia, The Free Encyclopedia. Retrieved 16:18, November 28, 2018,
           from https://en.wikipedia.org/w/index.php?title=Vectorization_(mathematics)&oldid=859035294
