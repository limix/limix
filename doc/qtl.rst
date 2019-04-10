************************
Quantitative trait locus
************************

Introduction
============

Every genetic model considered here is an instance of **generalized linear mixed model**
(GLMM).
It consists in four main components [St16]_:

- A linear predictor, 𝐳 = M𝛃 + 𝚇𝐮.
- The distribution of the random effects, 𝐮 ∼ 𝓝(𝟎, Σ).
- The residual distribution, yᵢ | 𝐮.
- The link function, g(𝜇ᵢ) = zᵢ.

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
of the identity link function, g(𝜇ᵢ) = 𝜇ᵢ, and of normally distributed residuals, yᵢ |
𝐮 ∼ 𝓝(𝜇ᵢ, 𝜎ᵢ²) [Mc11]_. It is more commonly described by the equation

.. math::
    :label: lmm

    𝐲 = 𝙼𝛃 + 𝚇𝐮 + 𝛆,

for which 𝜀ᵢ∼𝓝(0, 𝜎ᵢ²).  The random variables 𝐮 and 𝛆 are independent from each
other as well as 𝜀ᵢ and 𝜀ⱼ for i≠j.  Defining 𝐯 = 𝚇𝐮 leads to:

.. math::

    𝐯 ∼ 𝓝(𝟎, 𝚇Σ𝚇ᵀ).

There is another even simpler instance of GLMM that is also used in genetic analysis:
a **linear model** (LM) is merely a LMM without the random effects:

.. math::

    𝐲 = 𝙼𝛃 + 𝛆.

The above models are used to establish a statistical tests to find significant
association between genetic loci and phenotype. For that, their parameters have to be
estimated.

As an example, let us define two parameters that will describe the overall variances of
the random effects and of the residual effects:

.. math::

    Σ = 𝓋₀𝙸₀ ~~\text{and}~~ 𝜎ᵢ² = 𝓋₁.

If we assume a LMM, this example of model can be described by Eq. :eq:`lmm` for which

.. math::

    𝐯∼𝓝(𝟎, 𝓋₀𝚇𝚇ᵀ) ~~\text{and}~~ 𝛆∼𝓝(𝟎, 𝓋₁𝙸₁).

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
We will make use of the LRT approach in the next sections for flag significant genetic
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
    >>> print(r) # doctest: +FLOAT_CMP

The variable ``r`` is instance of the class :class:`limix.qtl.ScanResult` and stores all
the results of the analysis.  Printing it as we did above it will show a summary of the
results.

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

    >>> from limix.stats import linear_kinship
    >>> from numpy import zeros
    >>>
    >>> # Whole genotype of each sample.
    >>> X = random.randn(n, 50)
    >>> # Estimate a kinship relationship between samples.
    >>> K = linear_kinship(X, verbose=False)
    >>> # Update the phenotype
    >>> y += random.multivariate_normal(zeros(n), K)
    >>>
    >>> r = scan(X, y, "normal", K, 𝙼=M, verbose=False)
    >>> print(r) # doctest: +FLOAT_CMP
    Null model
    ----------
    <BLANKLINE>
    𝐲 ~ 𝓝(𝙼𝜶, 1.5632⋅𝙺 + 3.2301⋅𝙸)
    𝙼     = ['offset' 'age']
    𝜶     = [-1.88025701  0.19028836]
    se(𝜶) = [0.327493   0.01222069]
    lml   = -215.9781119592618
    <BLANKLINE>
    Alt model
    ---------
    <BLANKLINE>
    𝐲 ~ 𝓝(𝙼𝜶 + 𝙶𝞫, 1.5632⋅𝙺 + 3.2301⋅𝙸)
    min(pv)  = 0.014200670407257475
    max(lml) = -212.97159992879634
    <BLANKLINE>

Non-normal trait association
============================

If the residuals of the phenotype does not follow a Normal distribution, then we might
consider performing the analysis using a **generalized linear mixed model**. Let us
consider Poisson distributed residuals:

.. math::

    yᵢ | 𝐳 ∼ 𝙿𝚘𝚒𝚜𝚜𝚘𝚗(g(𝜇ᵢ)=zᵢ),

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
    >>> print(r) # doctest: +FLOAT_CMP
    Null model
    ----------
    <BLANKLINE>
    𝐳 ~ 𝓝(𝙼𝜶, 0.1130⋅𝙺 + 0.1399⋅𝙸) for yᵢ ~ Poisson(λᵢ=g(zᵢ)) and g(x)=eˣ
    𝙼     = ['offset' 'age']
    𝜶     = [-1.41641664  0.05496354]
    se(𝜶) = [0.2020572 0.0060997]
    lml   = -151.15802807711944
    <BLANKLINE>
    Alt model
    ---------
    <BLANKLINE>
    𝐳 ~ 𝓝(𝙼𝜶 + 𝙶𝞫, 0.1130⋅𝙺 + 0.1399⋅𝙸) for yᵢ ~ Poisson(λᵢ=g(zᵢ)) and g(x)=eˣ
    min(pv)  = 0.004370366590054564
    max(lml) = -147.09645533116503
    <BLANKLINE>

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

.. doctest::

    >>> from numpy import concatenate, newaxis
    >>> from limix.qtl import scan
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