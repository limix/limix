************************
Quantitative trait locus
************************


Introduction
^^^^^^^^^^^^

Every genetic model considered here is an instance of **generalised linear mixed model**
(GLMM).
It consists in four main components [St16]_:

    - A linear predictor, 𝐳 = M𝛃 + X𝐮.
    - The distribution of the random effects, 𝐮 ∼ 𝓝(𝟎, Σ).
    - The distribution of the outcome conditioned on the random effects (also known as
      the residual distribution), yᵢ | 𝐮.
    - The link function, g(𝜇ᵢ) = zᵢ.

The term 𝜇ᵢ represents the mean of yᵢ conditioned on 𝐮::

    𝜇ᵢ = \mathbb E[yᵢ | 𝐮].

The role of the link function is to scale the domain of zᵢ, which ranges from
-∞ to +∞, to the residual distribution parameter 𝜇ᵢ.
For example, the mean of a Bernoulli distribution is bounded within [0, 1], and
therefore requires a link function to translate values of zᵢ into values of
𝜇ᵢ.

The distribution of the outcome, conditioned on the random effects, has to be one from
the exponential family [Ef18]_ having mean 𝜇ᵢ::

    yᵢ | 𝐮 ∼ \text{ExpFam}(𝜇ᵢ).

A notable instance of the above model is the **linear mixed model** (LMM).
It consists of the identity link function, g(𝜇ᵢ) = 𝜇ᵢ, and of normally
distributed residuals, yᵢ | 𝐮 ∼ 𝓝(𝜇ᵢ, 𝜎ᵢ²)
[Mc11]_.
It is more commonly described by the equation ::

    𝐲 = M𝛃 + X𝐮 + 𝛆,

for which 𝜀ᵢ∼𝓝(0, 𝜎ᵢ²).  The random variables
𝐮 and 𝛆 are independent from each other as
well as 𝜀ᵢ and 𝜀ⱼ for i≠j.  Defining
𝐯 = X𝐮 leads to ::

    𝐯 ∼ 𝓝(𝟎, XΣXᵀ).

There is another even simpler instance of GLMM that is also used in genetic analysis:
a **linear model** (LM) is merely a LMM without the random effects::

    𝐲 = M𝛃 + 𝛆.

The above models are used to establish a statiscal tests to find significant association
between genetic loci and phenotype.
For that, their parameters have to be estimated.

As an example, let us define two parameters that will describe the overall variances of
the random effects
and of the residual effects::

    Σ = v₀I₀ and 𝜎ᵢ² = v₁.

If we assume a LMM, this example of model can be described by Eq. :eq:`eq_lmm` for which
::

    𝐯∼𝓝(𝟎, v₀XXᵀ) and 𝛆∼𝓝(𝟎, v₁I₁).

Therefore we have a model with three parameters: an array of effect sizes
𝛃
and variances v₀ and v₁.

Statistical test
^^^^^^^^^^^^^^^^

We use the **likelihood ratio test** (LRT) approach [LR18]_ to assess the significance
of the association
between genetic variants and the phenotype.
It is based on the ratio between the marginal likelihood of the null and alternative
models:

.. math::

    \mathcal H₀: 𝛉₀\\
    \mathcal H₁: 𝛉₁

where 𝛉₀ is fit under the null model and
𝛉₁ is fit under the alternative model.
The parameter inference is done via the maximum likelihood estimation (MLE) approach
[ML18]_::

    \boldsymbol{\hat{\theta}} = \underset{𝛉}{\mathrm{argmax}}~~
        p(𝐲 | M, X; 𝛉).

The likelihood ratio is then equal to ::

    \frac{p(𝐲| M, X; \boldsymbol{\hat{\theta₀}})}
        {p(𝐲| M, X; \boldsymbol{\hat{\theta₁}})}.

which will define the p-value of that comparison.

Single-trait association
^^^^^^^^^^^^^^^^^^^^^^^^

We first consider that the observed phenotype is described by additive effects from
covariates and genetic components, and any deviation from that is captured by the
assumed residual distribution and/or an over-dispersion component.  Let :math:`\mathbf
M` be a matrix of covariates and let \mathbf G be a matrix of genetic variants
that we suspect might have some effect on the phenotype.  Therefore, we have the linear
model::

    𝐲 = \underbrace{M\boldsymbol\alpha}_{\text{covariates}}+
    \underbrace{\mathbf G𝛃}_{\text{genetics}}+
    \underbrace{𝛆}_{\text{noise}},\\
    \text{where}~~𝛆∼𝓝(𝟎, v₁I),~~~~~~

and we wish to compare the following hypotheses::

    \mathcal H₀: 𝛃 = 0\\
    \mathcal H₁: 𝛃 ≠ 0

Note that the parameters of the above model are the covariate effect sizes,
\boldsymbol\alpha, the effect sizes of a set of genetic variants,
𝛃, and the variance v₁ of the noise variable.  Under the
null hypothesis, we set 𝛃=𝟎 and fit the rest of the
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

We now apply the function :func:`limix.qtl.scan` to our dataset

.. doctest::

    >>> from limix.qtl import scan
    >>>
    >>> r = scan(snps, y, 'normal', M=M, verbose=False)
    >>> print(r)
    Null model
    ----------
    <BLANKLINE>
      𝐲 ~ 𝓝(M𝜶, 0.32*K + 0.00*I)
      M = ['offset' 'age']
      𝜶 = [-0.81858684  0.02015968]
      Log marg. lik.: -21.218829574364268
      Number of models: 1
    <BLANKLINE>
    Alt model
    ---------
    <BLANKLINE>
      𝐲 ~ 𝓝(M𝜶 + Gᵢ, 0.32*K + 0.00*I)
      Min. p-value: 0.02219182245364262
      First perc. p-value: 0.0262094622393102
      Max. log marg. lik.: -18.60348830672571
      99th perc. log marg. lik.: -18.651776372344084
      Number of models: 4
    <BLANKLINE>

The variable ``r`` is instance of the class :class:`limix.qtl.QTLResult` and stores all
the results of the analysis.  Printing it as we did above it will show a summary of the
results.

Suppose we also have access to the whole genotype of our samples, X, and
we want to use them to account for population structure and cryptic relatedness in our
data (CITE).  Since the number of genetic variants in X is commonly
larger than the number of samples, and because we are not acctually interested in their
effect sizes, we will include it in our model as a random component.  We now have a
**linear mixed model**:

.. math::

    𝐲 = \underbrace{M\boldsymbol\alpha}_{\text{covariates}}+
    \underbrace{\mathbf G𝛃}_{\text{genetics}}+
    \underbrace{X𝐮}_{\text{pop. struct.}}+
    \underbrace{𝛆}_{\text{noise}},\\
    \text{where}~~
        𝐮∼𝓝(𝟎, v₀I₀) ~~\text{and}
    ~~𝛆∼𝓝(𝟎, v₁I₁).

It is important to note that 𝐯=X𝐮 can be equivalenty
described by a multivariate Normal distribution with a covariance proportional to
\mathbf K = XXᵀ::

    𝐯 ∼ 𝓝(𝟎, v₀\mathbf K).

We perform the analysis again now using also the covariance \mathbf K by calling
the function :func:`limix.qtl.scan`.

.. doctest::

    >>> from limix.stats import linear_kinship
    >>>
    >>> # Whole genotype of each sample.
    >>> X = random.randn(n, 50)
    >>> # Estimate a kinship relationship betweem samples.
    >>> K = linear_kinship(X, verbose=False)
    >>>
    >>> result = scan(X, y, 'normal', K, M=M, verbose=False)
    >>> print(result.stats.head()) # doctest: +FLOAT_CMP
          null lml   alt lml   pvalue  dof
    test
    0    -21.21883 -21.15531  0.72152    1
    1    -21.21883 -21.09391  0.61718    1
    2    -21.21883 -20.92358  0.44223    1
    3    -21.21883 -21.21649  0.94542    1
    4    -21.21883 -20.87087  0.40416    1
    >>> print(result.alt_effsizes.head()) # doctest: +FLOAT_CMP
       test candidate  effsize  effsize se
    0     0         0  0.04675     0.13116
    1     1         1 -0.05855     0.11713
    2     2         2 -0.09668     0.12582
    3     3         3  0.00746     0.10899
    4     4         4  0.12734     0.15264
    >>> print(result) # doctest: +FLOAT_CMP
    Null model
    ----------
    <BLANKLINE>
      𝐲 ~ 𝓝(M𝜶, 0.32*K + 0.00*I)
      M = ['offset' 'age']
      𝜶 = [-0.81858684  0.02015968]
      Log marg. lik.: -21.21882957624215
      Number of models: 1
    <BLANKLINE>
    Alt model
    ---------
    <BLANKLINE>
      𝐲 ~ 𝓝(M𝜶 + Gᵢ, 0.32*K + 0.00*I)
      Min. p-value: 0.01042644226036883
      First perc. p-value: 0.016787533334797423
      Max. log marg. lik.: -17.93855702329621
      99th perc. log marg. lik.: -18.28709258817481
      Number of models: 50

Generalised phenotype
~~~~~~~~~~~~~~~~~~~~~

If the residuals of the phenotype does not follow a Normal distribution, then we might
consider perform analysis using a **generalised linear mixed model**.  Let us consider
Poisson distributed residuals::

    yᵢ | 𝐳 ∼ \text{Bernoulli}(g(𝜇ᵢ)=zᵢ).

In the latter case, the 𝛆 can be used to describe the
dispersion between samples not fully captured by the residual distribution.

The following example applies :func:`limix.qtl.scan` to perform five likelihood ratio
tests for association with an outcome vector ``y`` having residual errors that follow a
Poisson distribution.  The matrix ``G`` defines both the five alternative hypotheses
(the first five columns) and the covariance matrix (the remaining columns).

.. doctest::

    >>> from numpy import exp, sqrt
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
    >>> K = linear_kinship(G[:, 5:], verbose=False)
    >>> result = scan(candidates, y, 'poisson', K, verbose=False)
    >>>
    >>> print(result.stats.head()) # doctest: +FLOAT_CMP
          null lml   alt lml   pvalue  dof
    test
    0    -34.64566 -33.88180  0.21645    1
    1    -34.64566 -34.35004  0.44194    1
    2    -34.64566 -34.40067  0.48394    1
    3    -34.64566 -33.93787  0.23413    1
    4    -34.64566 -34.56898  0.69534    1
    >>> print(result.alt_effsizes.head()) # doctest: +FLOAT_CMP
       test candidate  effsize  effsize se
    0     0         0  1.62727     1.31655
    1     1         1 -1.02366     1.33129
    2     2         2 -1.23573     1.76537
    3     3         3  1.97540     1.66030
    4     4         4 -0.53729     1.37198
    >>> print(result) # doctest: +FLOAT_CMP
    Null model
    ----------
    <BLANKLINE>
      𝐳 ~ 𝓝(M𝜶, 0.00*K + 0.03*I)
      yᵢ ~ Poisson(λᵢ=g(zᵢ)), where g(x)=eˣ
      M = ['offset']
      𝜶 = [-0.0141227]
      Log marg. lik.: -34.645664448446965
      Number of models: 1
    <BLANKLINE>
    Alt model
    ---------
    <BLANKLINE>
      𝐳 ~ 𝓝(M𝜶 + Gᵢ, 0.00*K + 0.03*I)
      yᵢ ~ Poisson(λᵢ=g(zᵢ)), where g(x)=eˣ
      Min. p-value: 0.21645253947712215
      First perc. p-value: 0.2171596825117883
      Max. log marg. lik.: -33.88179641668344
      99th perc. log marg. lik.: -33.88403939629015
      Number of models: 5

Single-trait with interaction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following linear mixed model is considered::

    \mathbf{y} =
    \underbrace{M𝛃}_
            {\substack{\text{fixed effects}\\ \text{without interaction}}}+
    \underbrace{(\mathbf G\odot\mathbf E₀)𝛃₀}_{\mathrm G\times\mathrm E₀} +
    \underbrace{\mathbf G\odot\mathbf E₁𝛃₁}_{\mathrm G\times\mathrm E₁} +
    \underbrace{X\mathbf{u}}_{\text{random effects}}+
    \underbrace{\boldsymbol{𝜀}}_{\text{residual}}.

The **GxE** terms are also fixed effects but encoding the interations between genetic
variants and environmental covariates defined by the user.

.. doctest::

    >>> from numpy import concatenate, newaxis
    >>> from limix.qtl import stᵢscan
    >>>
    >>> # generate interacting variables (environment)
    >>> random = RandomState(1)
    >>> E = random.randn(y.shape[0], 1)
    >>>
    >>> # add additive environment as covariate
    >>> ME = concatenate([M, E], axis=1)
    >>>
    >>> snps = random.randn(n, 100)
    >>>

# interaction test
res = stᵢscan(snps, y[:, newaxis], M=ME, E1=E, verbose=False)
print(res.head())  # doctest: +FLOAT_CMP
       pv1      pv0       pv    beta0  beta0_ste     lrt1     lrt0      lrt
0  0.14584  0.06186  0.54644  0.36731    0.19671  3.85044  3.48671  0.36373
1  0.81134  0.52514  0.90466  0.13994    0.22022  0.41813  0.40378  0.01435
2  0.74902  0.45885  0.86414  0.17079    0.23056  0.57798  0.54871  0.02928
3  0.77165  0.79650  0.50141 -0.05915    0.22937  0.51845  0.06650  0.45195
4  0.81176  0.64857  0.64725  0.09675    0.21229  0.41709  0.20770  0.20939


The process method returns three sets of P values: (i) ``pv0`` are association test P
values (\boldsymbol{\alpha}≠{0} when \boldsymbol{\beta}={0}), (ii)
``pv1`` are association + interaction P values (:math:`\left[\boldsymbol{\beta},
\boldsymbol{\alpha}\right]≠{0}`) and (iii) ``pv`` are interaction P values
(\boldsymbol{\alpha}≠{0}).  The effect sizes of the association test are also
returned.

If ``E0`` is not specified, a column-vector of ones is considered.  In this case the
\mathbf G\odot\mathbf E₀ term reduces to an additive genetic effect, and thus
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


# interaction test
r = stᵢscan(snps, y[:, newaxis], M=ME, E1=E1, E0=E0, verbose=False)
print(r.head())  # doctest: +FLOAT_CMP
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
    \underbrace{\mathbf{M}𝛃}_{\text{covariates}}+
    \underbrace{\mathbf{x}\odot\boldsymbol\gamma}_{\text{genetics}}+
    \underbrace{\mathbf E𝐮}_{\text{random effects}}+
    \underbrace{𝛆}_{\text{noise}},

where

.. math::
    \boldsymbol\gamma∼𝓝(𝟎,
    𝜎²_g(\underbrace{(1-\rho)\mathbf 1}_{\text{persistent}}
        + \underbrace{\rho\mathbf E\mathbf Eᵀ}_{\text{GxE}}),\\
    𝐮∼𝓝(𝟎, v₀I),
    ~~\text{and}~~
    𝛆∼𝓝(𝟎, v₁I).

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
    \underbrace{(\mathbf A_c \otimes M) \text{vec}(\mathbf B_c)}_{\text{covariates}}+
    \underbrace{(\mathbf A_g \otimes \mathbf G) \text{vec}(\mathbf B_g)}_{\text{genetics}}+
    \underbrace{\text{vec}(\mathbf U)}_{\text{random effect}}+
    \underbrace{\text{vec}(\boldsymbol\Psi)}_{\text{noise}}

is equivalent to Eq. :eq:`eq_lmm` but structured in a different way.
The columns of \mathbf Y correspond to the different traits being
considered.
The columns of \mathbf Y are stacked over each other and is denoted by
\text{vec}(\mathbf Y).
This is a linear transformation called vectorization [Ve19]_, and helps us
describe the model in a more concise manner.

The matrices \mathbf A_c and \mathbf A_g are design matrices for
the covariates and genetic variants, respectively.
The random effect component is defined by

.. math::

    \text{vec}(\mathbf U)∼𝓝(𝟎, \mathbf C₀\otimes\mathbf K)

and the residuals by

.. math::

    \text{vec}(\boldsymbol\Psi)∼𝓝(𝟎, \mathbf C₁\otimesI_n).

As before, M is the covariates matrix and \mathbf G is the
matrix of genetic variants.
The matrices \mathbf C₀ and \mathbf C₁ are two matrix-parameters and,
us such, are fitted during the likelihood maximisation.

Any-effect association test
~~~~~~~~~~~~~~~~~~~~~~~~~~~

An any-effect association test corresponds to testing
𝛃≠𝟎 with \mathbf A_g = I.

.. doctest::

    >>> from limix.qtl import scan
    >>> from numpy import eye
    >>>
    >>> p = 4
    >>> Y = random.randn(y.shape[0], p)
    >>>
    >>> Asnps = eye(p)
    >>> r = scan(G, Y, M=M, Asnps=Asnps, K=K, verbose=False)
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
    >>> r = scan(G, Y, K=K, Ac=None, Asnps=Asnps, Asnps0=Asnps0, verbose=False)
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
