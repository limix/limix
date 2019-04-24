************************
Quantitative trait locus
************************

Introduction
============

Every genetic model considered here is an instance of **generalized linear mixed model**
(GLMM).
It consists in four main components [St16]_:

- A linear predictor, 𝐳 = M𝛂 + 𝚇𝐮.
- The distribution of the random effects, 𝐮 ∼ 𝓝(𝟎, Σ).
- The residual distribution, yᵢ | 𝐮.
- The link function, 𝜇ᵢ = g(zᵢ).

The term 𝜇ᵢ represents the mean of yᵢ conditioned on 𝐮:

.. math::

    𝜇ᵢ = 𝙴[yᵢ|𝐮].

The role of the link function is to scale the domain of zᵢ, which ranges from -∞ to +∞,
to the residual distribution parameter 𝜇ᵢ. For example, the mean of a Bernoulli
distribution is bounded within [0, 1], and therefore requires a link function to
translate values of zᵢ into values of
𝜇ᵢ.

The distribution of the outcome, conditioned on the random effects, has to be one from
the exponential family [Ef18]_ having mean 𝜇ᵢ:

.. math::

    yᵢ|𝐮 ∼ 𝙴𝚡𝚙𝙵𝚊𝚖(𝜇ᵢ).

A notable instance of the above model is the **linear mixed model** (LMM). It consists
of the identity link function, 𝜇ᵢ = g(𝜇ᵢ), and of normally distributed residuals, yᵢ |
𝐮 ∼ 𝓝(𝜇ᵢ, 𝜎ᵢ²) [Mc11]_. It is more commonly described by the equation

.. math::
    :label: lmm

    𝐲 = 𝙼𝛂 + 𝚇𝐮 + 𝛆,

for which 𝜀ᵢ∼𝓝(0, 𝜎ᵢ²).  The random variables 𝐮 and 𝛆 are independent from each
other as well as 𝜀ᵢ and 𝜀ⱼ for i≠j.  Defining 𝐯 = 𝚇𝐮 leads to:

.. math::

    𝐯 ∼ 𝓝(𝟎, 𝚇Σ𝚇ᵀ).

There is another even simpler instance of GLMM that is also used in genetic analysis:
a **linear model** (LM) is merely a LMM without the random effects:

.. math::

    𝐲 = 𝙼𝛂 + 𝛆.

The above models are used to establish a statistical tests to find significant
association between genetic loci and phenotype. For that, their parameters have to be
estimated.

As an example, let us define two parameters that will describe the overall variances of
the random effects and of the residual effects:

.. math::

    Σ = 𝓋₀𝙸₀ ~~\text{and}~~ 𝜎ᵢ² = 𝓋₁.

If we assume a LMM, this example of model can be described by Eq. :eq:`lmm` for which

.. math::

    𝐮 ∼ 𝓝(𝟎, 𝓋₀𝙸₀) ~~\text{and}~~ 𝛆 ∼ 𝓝(𝟎, 𝓋₁𝙸₁).

Equivalently, we have

.. math::

    𝐲 = 𝙼𝛂 + 𝐯 + 𝛆,

for which

.. math::

    𝐯 ∼ 𝓝(𝟎, 𝓋₀𝚇𝚇ᵀ) ~~\text{and}~~ 𝛆 ∼ 𝓝(𝟎, 𝓋₁𝙸₁).

Therefore we have a model with three parameters: an array of effect sizes 𝛃 and
variances 𝓋₀ and 𝓋₁. If 𝚇 contains the normalized SNP genotypes of the samples, 𝚇𝚇ᵀ is
an estimation of the genetic relationship between the samples [Wa17]_.

Statistical test
================

We use the **likelihood ratio test** (LRT) approach [LR18]_ to assess the significance
of the association
between genetic variants and the phenotype.
It is based on the ratio between the marginal likelihood of the null 𝓗₀ and alternative
𝓗₁ models, for which the simpler model 𝓗₀ is defined by constraint one or more
parameters if the alternative model 𝓗₁.

The parameter inference is done via the maximum likelihood estimation (MLE) approach
[ML18]_, for which the marginal likelihood p(𝐲 | 𝙼, 𝚇; 𝛉) is maximized over the
parameters set 𝛉.
Let 𝛉₀ and 𝛉₁ be the optimal parameters set under the null and alternative models.
The likelihood ratio statistics is give by

.. math::

    -2 \log(p(𝐲| 𝙼, 𝚇; 𝛉₀) / p(𝐲| 𝙼, 𝚇; 𝛉₁)),

which asymptotically follows a χ² distribution [Wh14]_.
We will make use of the LRT approach in the next sections to flag significant genetic
associations.

Single-trait association
========================

We first consider that the observed phenotype is described by additive effects from
covariates and genetic components. Any deviation from that is assumed to be captured by
the residual distribution. Let 𝙼 be a matrix of covariates and let 𝙶 be a matrix of
genetic variants that we suspect might have some effect on the phenotype. Therefore, we
have the linear model:

.. math::

    𝐲 = \underbrace{𝙼𝛂}_{\text{covariates}}+
        \underbrace{𝙶𝛃}_{\text{genetics}}+
        \underbrace{𝛆}_{\text{noise}},\\
        \text{where}~~𝛆∼𝓝(𝟎, 𝓋₁𝙸),~~~~~~

and we wish to compare the following hypotheses:

.. math::

    𝓗₀: 𝛃 = 𝟎\\
    𝓗₁: 𝛃 ≠ 𝟎

Note that the parameters of the above model are the covariate effect sizes, 𝛂, the
effect sizes of a set of genetic variants, 𝛃, and the variance 𝓋₁ of the noise
variable.  Under the null hypothesis, we set 𝛃=𝟎 and fit the rest of the parameters.
Under the alternative hypothesis, we learn all the parameters. At the end, we compare
the marginal likelihoods via the likelihood ratio test.

Let us first generate a random data set having a phenotype, covariates, and a set of
genetic candidates.

.. doctest::

    >>> from numpy import ones, stack
    >>> from numpy.random import RandomState
    >>> from pandas import DataFrame
    >>>
    >>> random = RandomState(2)
    >>>
    >>> # sample size
    >>> n = 100
    >>>
    >>> # covariates
    >>> offset = ones(n) * random.randn()
    >>> age = random.randint(16, 75, n)
    >>> M = stack((offset, age), axis=1)
    >>> M = DataFrame(stack([offset, age], axis=1), columns=["offset", "age"])
    >>> M["sample"] = [f"sample{i}" for i in range(n)]
    >>> M = M.set_index("sample")
    >>> print(M.head())
              offset      age
    sample
    sample0 -0.41676 38.00000
    sample1 -0.41676 59.00000
    sample2 -0.41676 34.00000
    sample3 -0.41676 27.00000
    sample4 -0.41676 56.00000
    >>> # genetic variants
    >>> G = random.randn(n, 4)
    >>>
    >>> # sampling the phenotype
    >>> alpha = random.randn(2)
    >>> beta = random.randn(4)
    >>> eps = random.randn(n)
    >>> y = M @ alpha + G @ beta + eps

We now apply the function :func:`limix.qtl.scan` to our data set

.. doctest::

    >>> from limix.qtl import scan
    >>>
    >>> r = scan(G, y, "normal", M=M, verbose=False)
    >>> print(r)
    Hypothesis 0
    ============
    <BLANKLINE>
    𝐲 ~ 𝓝(𝙼𝜶, 3.462⋅𝙸)
    <BLANKLINE>
    M     = ['offset' 'age']
    𝜶     = [2.10096551 0.19582931]
    se(𝜶) = [1.25826998 0.01068367]
    lml   = -203.98750767964498
    <BLANKLINE>
    Hypothesis 2
    ============
    <BLANKLINE>
    𝐲 ~ 𝓝(𝙼𝜶 + G𝛃, s(3.462⋅𝙸))
    <BLANKLINE>
              lml       cov. effsizes   cand. effsizes
    --------------------------------------------------
    mean   -1.951e+02       9.915e-01       -6.198e-01
    std     9.227e+00       9.342e-01        3.974e-01
    min    -2.031e+02       1.844e-01       -1.025e+00
    25%    -2.026e+02       1.959e-01       -9.275e-01
    50%    -1.967e+02       5.965e-01       -6.047e-01
    75%    -1.893e+02       1.831e+00       -2.970e-01
    max    -1.841e+02       2.312e+00       -2.448e-01
    <BLANKLINE>
    Likelihood-ratio test p-values
    ==============================
    <BLANKLINE>
           𝓗₀ vs 𝓗₂
    ----------------
    mean   6.514e-02
    std    8.856e-02
    min    2.804e-10
    25%    2.606e-07
    50%    3.651e-02
    75%    1.016e-01
    max    1.875e-01

Suppose we also have access to the whole genotype of our samples, 𝚇, and we want to use
them to account for population structure and cryptic relatedness in our data [Ho13]_.
Since the number of genetic variants in 𝚇 is commonly larger than the number of
samples, and because we are not actually interested in their effect sizes, we will
include it in our model as a random component. We now have a **linear mixed model**:

.. math::

    𝐲 = \underbrace{𝙼𝛂}_{\text{covariates}}+
        \underbrace{𝙶𝛃}_{\text{genetics}}+
        \underbrace{𝚇𝐮}_{\text{pop. struct.}}+
        \underbrace{𝛆}_{\text{noise}},\\
        \text{where}~~
            𝐮∼𝓝(𝟎, 𝓋₀𝙸₀) ~~\text{and}
            ~~𝛆∼𝓝(𝟎, 𝓋₁𝙸₁).

It is important to note that 𝐯=𝚇𝐮 can be equivalently described by a multivariate
Normal distribution with a covariance proportional to 𝙺 = 𝚇𝚇ᵀ:

.. math::

    𝐯 ∼ 𝓝(𝟎, 𝓋₀𝙺).

We make use of the function :func:`limix.stats.linear_kinship` to define the covariance
matrix 𝙺, and call :func:`limix.qtl.scan` to perform the analysis.

.. doctest::

    >>> from limix.stats import linear_kinship, multivariate_normal
    >>> from numpy import zeros, eye
    >>>
    >>> # Whole genotype of each sample.
    >>> X = random.randn(n, 50)
    >>> # Estimate a kinship relationship between samples.
    >>> K = linear_kinship(X, verbose=False) + 1e-9 * eye(n)
    >>> # Update the phenotype
    >>> y += multivariate_normal(random, zeros(n), K)
    >>>
    >>> r = scan(X, y, "normal", K, 𝙼=M, verbose=False)
    >>> print(r)
    Hypothesis 0
    ============
    <BLANKLINE>
    𝐲 ~ 𝓝(𝙼𝜶, 1.436⋅𝙺 + 2.934⋅𝙸)
    <BLANKLINE>
    M     = ['offset' 'age']
    𝜶     = [1.95338293 0.19448903]
    se(𝜶) = [1.25455536 0.0107647 ]
    lml   = -211.3819625136375
    <BLANKLINE>
    Hypothesis 2
    ============
    <BLANKLINE>
    𝐲 ~ 𝓝(𝙼𝜶 + G𝛃, s(1.436⋅𝙺 + 2.934⋅𝙸))
    <BLANKLINE>
              lml       cov. effsizes   cand. effsizes
    --------------------------------------------------
    mean   -2.109e+02       1.069e+00        5.922e-02
    std     7.210e-01       8.819e-01        2.474e-01
    min    -2.114e+02       1.919e-01       -5.204e-01
    25%    -2.113e+02       1.944e-01       -1.102e-01
    50%    -2.111e+02       8.904e-01        4.920e-02
    75%    -2.109e+02       1.956e+00        2.366e-01
    max    -2.076e+02       2.215e+00        6.433e-01
    <BLANKLINE>
    Likelihood-ratio test p-values
    ==============================
    <BLANKLINE>
           𝓗₀ vs 𝓗₂
    ----------------
    mean   4.843e-01
    std    2.705e-01
    min    6.294e-03
    25%    3.137e-01
    50%    4.752e-01
    75%    6.953e-01
    max    9.929e-01

Non-normal trait association
============================

If the residuals of the phenotype does not follow a Normal distribution, then we might
consider performing the analysis using a **generalized linear mixed model**. Let us
consider Poisson distributed residuals:

.. math::

    yᵢ | 𝐳 ∼ 𝙿𝚘𝚒𝚜𝚜𝚘𝚗(𝜇ᵢ=g(zᵢ)),

where the latent phenotype is described by

.. math::

    𝐳 = 𝙼𝛃 + 𝚇𝐮 + 𝛆,

for

.. math::

    𝐮 ∼ 𝓝(𝟎, 𝓋₀𝙸₀) ~~\text{and}~~ 𝛆 ∼ 𝓝(𝟎, 𝓋₁𝙸₁).

Note that the term 𝛆 in the above model is not the residual variable, as it were in the
Eq. :eq:`lmm`.
The term 𝛆 is used to account for the so-called over-dispersion, i.e., when the residual
distribution is not sufficient to explain the variability of yᵢ.

.. doctest::

    >>> from numpy import exp
    >>>
    >>> z = (y - y.mean()) / y.std()
    >>> y = random.poisson(exp(z))
    >>>
    >>> r = scan(G, y, "poisson", K, M=M, verbose=False)
    >>> print(r)
    Hypothesis 0
    ============
    <BLANKLINE>
    𝐳 ~ 𝓝(𝙼𝜶, 0.154⋅𝙺 + 0.000⋅𝙸) for yᵢ ~ Poisson(λᵢ=g(zᵢ)) and g(x)=eˣ
    <BLANKLINE>
    M     = ['offset' 'age']
    𝜶     = [5.17511934 0.04665214]
    se(𝜶) = [0.85159296 0.00604329]
    lml   = -145.33385788740767
    <BLANKLINE>
    Hypothesis 2
    ============
    <BLANKLINE>
    𝐳 ~ 𝓝(𝙼𝜶 + G𝛃, s(0.154⋅𝙺 + 0.000⋅𝙸)) for yᵢ ~ Poisson(λᵢ=g(zᵢ)) and g(x)=eˣ
    <BLANKLINE>
              lml       cov. effsizes   cand. effsizes
    --------------------------------------------------
    mean   -1.440e+02       2.553e+00       -1.306e-01
    std     1.343e+00       2.682e+00        9.268e-02
    min    -1.453e+02       4.345e-02       -2.227e-01
    25%    -1.450e+02       4.635e-02       -2.018e-01
    50%    -1.439e+02       2.456e+00       -1.344e-01
    75%    -1.428e+02       5.054e+00       -6.321e-02
    max    -1.427e+02       5.202e+00       -3.085e-02
    <BLANKLINE>
    Likelihood-ratio test p-values
    ==============================
    <BLANKLINE>
           𝓗₀ vs 𝓗₂
    ----------------
    mean   2.830e-01
    std    3.213e-01
    min    2.274e-02
    25%    2.519e-02
    50%    2.113e-01
    75%    4.692e-01
    max    6.867e-01

Single-trait with interaction
=============================

The following linear mixed model is considered:

.. math::

    𝐲 = 𝙼𝛂 + (𝙶⊙𝙴₀)𝛃₀ + (𝙶⊙𝙴₁)𝛃₁ + 𝚇𝐮 + 𝛆,\\
    \text{where}~~ 𝐮∼𝓝(𝟎, 𝓋₀𝙸₀) ~~\text{and}~~ 𝛆∼𝓝(𝟎, 𝓋₁𝙸₁).

The operator ⊙ works as follows:

.. math::

    𝙰⊙𝙱 = [𝙰₀𝙱₀ ~~...~~ 𝙰₀𝙱ₙ ~~ 𝙰₁𝙱₀ ~~...~~ 𝙰₁𝙱ₙ ~~...~~ 𝙰ₘ𝙱ₙ]

Therefore, the terms 𝙶⊙𝙴₀ and 𝙶⊙𝙴₁ can be understood as interaction terms between
genetics, 𝙶, and environments, 𝙴₀ and 𝙴₁.

We define three hypotheses from the above linear mixed model:

.. math::

    𝓗₀: 𝛃₀=𝟎 ~~\text{and}~~ 𝛃₁=𝟎\\
    𝓗₁: 𝛃₀≠𝟎 ~~\text{and}~~ 𝛃₁=𝟎\\
    𝓗₂: 𝛃₀≠𝟎 ~~\text{and}~~ 𝛃₁≠𝟎

The hypothesis 𝓗₀ is for no-interaction, 𝓗₁ is for interaction with environments
encoded in 𝙴₀, and 𝓗₂ is for interaction with environments encoded in 𝙴₀ and 𝙴₁.
We perform three statistical tests:

- 𝓗₀ (null) vs 𝓗₁ (alternative)
- 𝓗₀ (null) vs 𝓗₂ (alternative)
- 𝓗₁ (null) vs 𝓗₂ (alternative)

Here is an example.

.. doctest::

    >>> from numpy import concatenate, newaxis
    >>> from limix.qtl import iscan
    >>>
    >>> # Generate interacting variables (environment)
    >>> E0 = random.randn(y.shape[0], 1)
    >>> E1 = random.randn(y.shape[0], 1)
    >>>
    >>> r = iscan(G, y, "normal", K, M, E0=E0, E1=E1, verbose=False)
    >>> print(r)
    Hypothesis 0
    ============
    <BLANKLINE>
    𝐲 ~ 𝓝(𝙼𝜶, 0.376⋅𝙺 + 2.077⋅𝙸)
    <BLANKLINE>
    M     = ['offset' 'age']
    𝜶     = [3.12608063 0.06042316]
    se(𝜶) = [1.01867609 0.00870181]
    lml   = -185.77488727691096
    <BLANKLINE>
    Hypothesis 1
    ============
    <BLANKLINE>
    𝐲 ~ 𝓝(𝙼𝜶 + (𝙶⊙𝙴₀)𝛃₀, s(0.376⋅𝙺 + 2.077⋅𝙸))
    <BLANKLINE>
              lml       cov. effsizes   cand. effsizes
    --------------------------------------------------
    mean   -1.856e+02       1.611e+00       -2.976e-03
    std     1.949e-01       1.658e+00        1.208e-01
    min    -1.858e+02       6.034e-02       -1.461e-01
    25%    -1.858e+02       6.058e-02       -4.769e-02
    50%    -1.856e+02       1.590e+00       -7.487e-03
    75%    -1.854e+02       3.137e+00        3.722e-02
    max    -1.854e+02       3.235e+00        1.492e-01
    <BLANKLINE>
    Hypothesis 1
    ============
    <BLANKLINE>
    𝐲 ~ 𝓝(𝙼𝜶 + (𝙶⊙𝙴₀)𝛃₀ + (𝙶⊙𝙴₁)𝛃₁, s(0.376⋅𝙺 + 2.077⋅𝙸))
    <BLANKLINE>
              lml       cov. effsizes   cand. effsizes
    --------------------------------------------------
    mean   -1.852e+02       1.612e+00        7.001e-03
    std     7.598e-01       1.659e+00        1.475e-01
    min    -1.857e+02       5.991e-02       -2.573e-01
    25%    -1.856e+02       6.096e-02       -4.135e-02
    50%    -1.855e+02       1.571e+00        3.611e-02
    75%    -1.851e+02       3.135e+00        7.660e-02
    max    -1.841e+02       3.241e+00        1.971e-01
    <BLANKLINE>
    Likelihood-ratio test p-values
    ==============================
    <BLANKLINE>
           𝓗₀ vs 𝓗₁    𝓗₀ vs 𝓗₂    𝓗₁ vs 𝓗₂
    ----------------------------------------
    mean   6.867e-01   6.501e-01   5.244e-01
    std    3.199e-01   3.350e-01   3.168e-01
    min    3.963e-01   1.795e-01   9.940e-02
    25%    4.185e-01   5.578e-01   3.784e-01
    50%    6.755e-01   7.277e-01   5.971e-01
    75%    9.436e-01   8.200e-01   7.431e-01
    max    9.995e-01   9.654e-01   8.042e-01


Multi-trait association
=======================

LMM can also be used to jointly model multiple traits.
Let n, c, and p be the number of samples, covariates, and traits, respectively.
The outcome variable 𝚈 is a n×p matrix distributed according to

..  math ::
    :label: mtlmm

    𝚟𝚎𝚌(𝚈) ∼ 𝓝((𝙰 ⊗ 𝙼) 𝚟𝚎𝚌(𝐀), 𝙲₀ ⊗ 𝚇𝚇ᵀ + 𝙲₁ ⊗ 𝙸).

𝙰 and 𝙼 are design matrices of dimensions p×p and n×c provided by the user,
where 𝙼 is the usual matrix of covariates commonly used in single-trait models.
𝐀 is a c×p matrix of fixed-effect sizes per trait.
𝚇 is a n×r matrix provided by the user and I is a n×n identity matrices.
𝙲₀ and 𝙲₁ are both symmetric matrices of dimensions p×p, for which 𝙲₁ is
guaranteed by our implementation to be of full rank.
The parameters of this model are the matrices 𝐀, 𝙲₀, and 𝙲₁.
𝚟𝚎𝚌(⋅) is a function that stacks the columns of the provided matrix into a vector
[Ve19]_.

Let 𝐲=𝚟𝚎𝚌(𝚈) and 𝛂=𝚟𝚎𝚌(𝐀).
We can extend the model in Eq. :eq:`mtlmm` to represent three different hypotheses:

..  math ::

    𝐲 ∼ 𝓝((𝙰 ⊗ 𝙼)𝛂 + (𝙰₀ ⊗ 𝙶)𝛃₀ + (𝙰₁ ⊗ 𝙶)𝛃₁, 𝙲₀ ⊗ 𝚇𝚇ᵀ + 𝙲₁ ⊗ 𝙸);

the hypotheses being

.. math::

    𝓗₀: 𝛃₀=𝟎 ~~\text{and}~~ 𝛃₁=𝟎\\
    𝓗₁: 𝛃₀≠𝟎 ~~\text{and}~~ 𝛃₁=𝟎\\
    𝓗₂: 𝛃₀≠𝟎 ~~\text{and}~~ 𝛃₁≠𝟎

as before.
Here is an example.

.. doctest::

    >>> from numpy import eye
    >>>
    >>> p = 2
    >>> Y = random.randn(n, p)
    >>> A = random.randn(p, p)
    >>> A = A @ A.T
    >>> A0 = ones((p, 1))
    >>> A1 = eye(p)
    >>>
    >>> r = scan(G, Y, K=K, M=M, A=A, A0=A0, A1=A1, verbose=False)
    >>> print(r)
    Hypothesis 0
    ============
    <BLANKLINE>
    𝐲 ~ 𝓝((A⊗𝙼)𝛂, C₀⊗𝙺 + C₁⊗𝙸)
    <BLANKLINE>
    traits   = ['0' '1']
    M        = ['offset' 'age']
    𝜶        = [-0.16350063 -0.00299804 -0.34519932 -0.00080396]
    se(𝜶)    = [11.30610534  0.09640495  5.36108799  0.04573768]
    diag(C₀) = [0.01405506 0.29149252]
    diag(C₁) = [0.81178933 0.85777798]
    lml      = -277.3341774005235
    <BLANKLINE>
    Hypothesis 1
    ============
    <BLANKLINE>
    𝐲 ~ 𝓝((A⊗𝙼)𝛂 + (A₀⊗G)𝛃₀, s(C₀⊗𝙺 + C₁⊗𝙸))
    <BLANKLINE>
              lml       cov. effsizes   cand. effsizes
    --------------------------------------------------
    mean   -2.763e+02      -1.243e-01       -2.842e-02
    std     1.329e+00       1.435e-01        1.120e-01
    min    -2.773e+02      -3.712e-01       -1.666e-01
    25%    -2.772e+02      -2.032e-01       -7.187e-02
    50%    -2.767e+02      -6.896e-02       -2.672e-02
    75%    -2.758e+02      -1.371e-03        1.673e-02
    max    -2.744e+02       1.388e-04        1.063e-01
    <BLANKLINE>
    Hypothesis 2
    ============
    <BLANKLINE>
    𝐲 ~ 𝓝((A⊗𝙼)𝛂 + (A₀⊗G)𝛃₀ + (A₁⊗G)𝛃₁, s(C₀⊗𝙺 + C₁⊗𝙸))
    <BLANKLINE>
              lml       cov. effsizes   cand. effsizes
    --------------------------------------------------
    mean   -2.761e+02      -1.245e-01       -1.209e-02
    std     1.404e+00       1.439e-01        5.823e-02
    min    -2.772e+02      -3.744e-01       -1.151e-01
    25%    -2.770e+02      -2.025e-01       -3.202e-02
    50%    -2.766e+02      -6.702e-02       -7.899e-03
    75%    -2.757e+02      -1.441e-03        2.127e-02
    max    -2.741e+02      -7.170e-04        7.371e-02
    <BLANKLINE>
    Likelihood-ratio test p-values
    ==============================
    <BLANKLINE>
           𝓗₀ vs 𝓗₁    𝓗₀ vs 𝓗₂    𝓗₁ vs 𝓗₂
    ----------------------------------------
    mean   3.973e-01   6.122e-01   8.438e-01
    std    3.851e-01   3.942e-01   1.327e-01
    min    1.597e-02   9.168e-02   7.251e-01
    25%    1.133e-01   4.159e-01   7.320e-01
    50%    3.626e-01   7.039e-01   8.370e-01
    75%    6.466e-01   9.002e-01   9.488e-01
    max    8.478e-01   9.493e-01   9.760e-01

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
.. [Wa17]  Wang, B., Sverdlov, S., & Thompson, E. (2017). Efficient estimation of
           realized kinship from single nucleotide polymorphism genotypes. Genetics,
           205(3), 1063-1078.
.. [Wh14]  White, H. (2014). Asymptotic theory for econometricians. Academic press.
.. [Ho13]  Hoffman, G. E. (2013). Correcting for population structure and kinship using
           the linear mixed model: theory and extensions. PloS one, 8(10), e75707.
