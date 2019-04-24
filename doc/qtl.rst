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
    >>> random = RandomState(1)
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
    sample0  1.62435 25.00000
    sample1  1.62435 27.00000
    sample2  1.62435 21.00000
    sample3  1.62435 31.00000
    sample4  1.62435 16.00000
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
    𝐲 ~ 𝓝(𝙼𝜶, 4.115⋅𝙸)
    <BLANKLINE>
    M     = ['offset' 'age']
    𝜶     = [-1.60130331  0.17922863]
    se(𝜶) = [0.33382518 0.01227417]
    lml   = -212.62741096350612
    <BLANKLINE>
    Hypothesis 2
    ============
    <BLANKLINE>
    𝐲 ~ 𝓝(𝙼𝜶 + G𝛃, s(4.115⋅𝙸))
    <BLANKLINE>
              lml       cov. effsizes   cand. effsizes
    --------------------------------------------------
    mean   -1.965e+02      -7.056e-01       -6.597e-01
    std     3.037e+01       9.470e-01        8.287e-01
    min    -2.125e+02      -1.648e+00       -1.891e+00
    25%    -2.118e+02      -1.585e+00       -7.375e-01
    50%    -2.112e+02      -6.789e-01       -3.283e-01
    75%    -1.959e+02       1.793e-01       -2.505e-01
    max    -1.509e+02       1.838e-01       -9.122e-02
    <BLANKLINE>
    Likelihood-ratio test p-values
    ==============================
    <BLANKLINE>
           𝓗₀ vs 𝓗₂
    ----------------
    mean   2.206e-01
    std    3.160e-01
    min    1.139e-28
    25%    4.727e-02
    50%    9.740e-02
    75%    2.707e-01
    max    6.876e-01

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
    𝐲 ~ 𝓝(𝙼𝜶, 1.809⋅𝙺 + 3.194⋅𝙸)
    <BLANKLINE>
    M     = ['offset' 'age']
    𝜶     = [-1.94732565  0.19294746]
    se(𝜶) = [0.32847297 0.01227062]
    lml   = -217.17874209561472
    <BLANKLINE>
    Hypothesis 2
    ============
    <BLANKLINE>
    𝐲 ~ 𝓝(𝙼𝜶 + G𝛃, s(1.809⋅𝙺 + 3.194⋅𝙸))
    <BLANKLINE>
              lml       cov. effsizes   cand. effsizes
    --------------------------------------------------
    mean   -2.167e+02      -8.781e-01        1.065e-02
    std     6.507e-01       1.077e+00        2.876e-01
    min    -2.172e+02      -2.047e+00       -6.198e-01
    25%    -2.172e+02      -1.949e+00       -1.781e-01
    50%    -2.170e+02      -8.392e-01        5.865e-03
    75%    -2.165e+02       1.930e-01        1.598e-01
    max    -2.150e+02       1.960e-01        6.612e-01
    <BLANKLINE>
    Likelihood-ratio test p-values
    ==============================
    <BLANKLINE>
           𝓗₀ vs 𝓗₂
    ----------------
    mean   5.185e-01
    std    3.223e-01
    min    3.889e-02
    25%    2.303e-01
    50%    5.590e-01
    75%    8.491e-01
    max    9.825e-01

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
    𝐳 ~ 𝓝(𝙼𝜶, 0.124⋅𝙺 + 0.096⋅𝙸) for yᵢ ~ Poisson(λᵢ=g(zᵢ)) and g(x)=eˣ
    <BLANKLINE>
    M     = ['offset' 'age']
    𝜶     = [-0.95485831  0.03996287]
    se(𝜶) = [0.16946147 0.00531426]
    lml   = -136.23239573615007
    <BLANKLINE>
    Hypothesis 2
    ============
    <BLANKLINE>
    𝐳 ~ 𝓝(𝙼𝜶 + G𝛃, s(0.124⋅𝙺 + 0.096⋅𝙸)) for yᵢ ~ Poisson(λᵢ=g(zᵢ)) and g(x)=eˣ
    <BLANKLINE>
              lml       cov. effsizes   cand. effsizes
    --------------------------------------------------
    mean   -1.349e+02      -4.575e-01       -8.332e-02
    std     2.145e+00       5.319e-01        1.631e-01
    min    -1.362e+02      -9.733e-01       -3.052e-01
    25%    -1.360e+02      -9.511e-01       -1.473e-01
    50%    -1.358e+02      -4.514e-01       -5.062e-02
    75%    -1.347e+02       3.996e-02        1.338e-02
    max    -1.317e+02       4.027e-02        7.322e-02
    <BLANKLINE>
    Likelihood-ratio test p-values
    ==============================
    <BLANKLINE>
           𝓗₀ vs 𝓗₂
    ----------------
    mean   4.206e-01
    std    3.920e-01
    min    2.559e-03
    25%    2.155e-01
    50%    3.713e-01
    75%    5.764e-01
    max    9.373e-01

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
    𝐲 ~ 𝓝(𝙼𝜶, 0.400⋅𝙺 + 2.242⋅𝙸)
    <BLANKLINE>
    M     = ['offset' 'age']
    𝜶     = [-0.62334332  0.06181413]
    se(𝜶) = [0.26023726 0.00964886]
    lml   = -189.46374093761898
    <BLANKLINE>
    Hypothesis 1
    ============
    <BLANKLINE>
    𝐲 ~ 𝓝(𝙼𝜶 + (𝙶⊙𝙴₀)𝛃₀, s(0.400⋅𝙺 + 2.242⋅𝙸))
    <BLANKLINE>
              lml       cov. effsizes   cand. effsizes
    --------------------------------------------------
    mean   -1.884e+02      -2.794e-01       -5.682e-02
    std     7.613e-01       3.649e-01        2.830e-01
    min    -1.895e+02      -6.345e-01       -3.100e-01
    25%    -1.888e+02      -6.195e-01       -2.423e-01
    50%    -1.882e+02      -2.722e-01       -1.231e-01
    75%    -1.879e+02       6.164e-02        6.247e-02
    max    -1.878e+02       6.243e-02        3.288e-01
    <BLANKLINE>
    Hypothesis 1
    ============
    <BLANKLINE>
    𝐲 ~ 𝓝(𝙼𝜶 + (𝙶⊙𝙴₀)𝛃₀ + (𝙶⊙𝙴₁)𝛃₁, s(0.400⋅𝙺 + 2.242⋅𝙸))
    <BLANKLINE>
              lml       cov. effsizes   cand. effsizes
    --------------------------------------------------
    mean   -1.880e+02      -2.790e-01       -4.654e-02
    std     6.498e-01       3.645e-01        2.302e-01
    min    -1.888e+02      -6.383e-01       -3.315e-01
    25%    -1.883e+02      -6.229e-01       -2.106e-01
    50%    -1.880e+02      -2.656e-01       -9.740e-02
    75%    -1.877e+02       6.165e-02        1.155e-01
    max    -1.872e+02       6.255e-02        3.352e-01
    <BLANKLINE>
    Likelihood-ratio test p-values
    ==============================
    <BLANKLINE>
           𝓗₀ vs 𝓗₁    𝓗₀ vs 𝓗₂    𝓗₁ vs 𝓗₂
    ----------------------------------------
    mean   3.035e-01   2.714e-01   3.901e-01
    std    3.935e-01   1.728e-01   1.688e-01
    min    6.780e-02   1.090e-01   2.482e-01
    25%    7.538e-02   1.682e-01   2.831e-01
    50%    1.286e-01   2.342e-01   3.423e-01
    75%    3.567e-01   3.374e-01   4.493e-01
    max    8.888e-01   5.084e-01   6.274e-01


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
    𝜶        = [-0.35658972  0.00603165 -0.27791642  0.00466769]
    se(𝜶)    = [0.39141908 0.01439179 0.31766476 0.01167998]
    diag(C₀) = [3.60751606e-12 2.59476377e-14]
    diag(C₁) = [0.89913313 1.14019828]
    lml      = -285.02921214513117
    <BLANKLINE>
    Hypothesis 1
    ============
    <BLANKLINE>
    𝐲 ~ 𝓝((A⊗𝙼)𝛂 + (A₀⊗G)𝛃₀, s(C₀⊗𝙺 + C₁⊗𝙸))
    <BLANKLINE>
              lml       cov. effsizes   cand. effsizes
    --------------------------------------------------
    mean   -2.845e+02      -1.543e-01        2.297e-02
    std     2.505e-01       1.684e-01        8.000e-02
    min    -2.849e+02      -3.846e-01       -9.417e-02
    25%    -2.846e+02      -3.066e-01        4.629e-03
    50%    -2.845e+02      -1.238e-01        5.502e-02
    75%    -2.844e+02       5.685e-03        7.336e-02
    max    -2.843e+02       7.323e-03        7.600e-02
    <BLANKLINE>
    Hypothesis 2
    ============
    <BLANKLINE>
    𝐲 ~ 𝓝((A⊗𝙼)𝛂 + (A₀⊗G)𝛃₀ + (A₁⊗G)𝛃₁, s(C₀⊗𝙺 + C₁⊗𝙸))
    <BLANKLINE>
              lml       cov. effsizes   cand. effsizes
    --------------------------------------------------
    mean   -2.838e+02      -1.542e-01        1.175e-02
    std     1.463e+00       1.684e-01        8.129e-02
    min    -2.849e+02      -3.841e-01       -1.290e-01
    25%    -2.844e+02      -3.080e-01       -3.439e-02
    50%    -2.843e+02      -1.215e-01        1.122e-02
    75%    -2.836e+02       5.721e-03        4.966e-02
    max    -2.816e+02       7.297e-03        1.898e-01
    <BLANKLINE>
    Likelihood-ratio test p-values
    ==============================
    <BLANKLINE>
           𝓗₀ vs 𝓗₁    𝓗₀ vs 𝓗₂    𝓗₁ vs 𝓗₂
    ----------------------------------------
    mean   3.587e-01   6.012e-01   6.956e-01
    std    1.632e-01   3.712e-01   4.387e-01
    min    2.279e-01   7.680e-02   5.880e-02
    25%    2.644e-01   5.317e-01   5.790e-01
    50%    3.062e-01   6.880e-01   8.620e-01
    75%    4.005e-01   7.576e-01   9.787e-01
    max    5.944e-01   9.521e-01   9.996e-01

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
