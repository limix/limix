************************
Quantitative trait locus
************************


Introduction
^^^^^^^^^^^^

Every genetic model considered here is an instance of **generalised linear mixed model** (GLMM).
It consists in four main components [St16]_:

    - A linear predictor, :math:`\mathbf z = \mathbf M\boldsymbol\beta + \mathbf X\mathbf u`.
    - The distribution of the random effects, :math:`\mathbf u \sim \mathcal N(\mathbf 0, \boldsymbol\Sigma)`.
    - The distribution of the outcome conditioned on the random effects, :math:`y_i | \mathbf u`. It is also
      sometimes called the residual distribution.
    - The link function, :math:`g(\mu_i) = z_i`.

The term :math:`\mu_i` represents the mean of :math:`y_i` conditioned on :math:`\mathbf u`:

.. math::

    \mu_i = \mathbb E[y_i | \mathbf u].

And the distribution of the outcome, conditioned on the random effects, is from the
exponential family [Ef18]_ having mean :math:`\mu_i`:

.. math::

    y_i | \mathbf u \sim \text{ExpFam}(\mu_i).

A notable instance of the above model is the **linear mixed model** (LMM).
It consist of the identity link function, :math:`g(\mu_i) = \mu_i`, and of normally distributed
residuals, :math:`y_i | \mathbf u \sim \mathcal N(\mu_i, \sigma_i^2)` [Mc11]_.
It is more commonly described by the equation

.. math::

    \mathbf y = \mathbf M\boldsymbol\beta + \mathbf X\mathbf u + \boldsymbol\epsilon,

for which :math:`\epsilon_i\sim\mathcal N(0, \sigma_i^2)`.
The random variables :math:`\mathbf u` and :math:`\boldsymbol\epsilon` are independent
from each other as well as :math:`\epsilon_i` and :math:`\epsilon_j` for :math:`i\neq j`.

A **linear model** (LM) is a LMM without the random component:

.. math::

    \mathbf y = \mathbf M\boldsymbol\beta + \boldsymbol\epsilon.

The above models are used to establish a statiscal test to find significant association between
genetic loci and phenotype.
Let
the null and alternative hypotheses:

.. math::

    \mathcal H_0:

Single-trait association
^^^^^^^^^^^^^^^^^^^^^^^^

We consider that the observed phenotype is described by additive effects from covariates
and genetic components, and any deviation from that is captured by the assumed residual
distribution and/or an over-dispersion component.
Additionally, the random component is often used to account for population structure and cryptic
relatedness.

As an example, let the matrix :math:`\mathbf M` consist of covariates (age, sex, etc.) and a few number
of genetic variants potentially associated with a given phenotype.
Let the matrix :math:`\mathbf X` consist of all genetic variants.
We assume that the phenotype under analysis is well described by

.. math::

    \mathbf z = \underbrace{\mathbf M\boldsymbol\beta}_{\text{fixed effects}}
        + \underbrace{\mathbf X\mathbf u}_{\text{random effects}}
        +  \underbrace{\boldsymbol\epsilon}_{
            \substack{\text{residual or}\\ \text{overdispersion}}
        }
        ,~\text{ for}\\
           \mathbf u\sim\mathcal N(\mathbf 0, v_0\mathrm I_0)
          ~~\text{and}~~
            \boldsymbol\epsilon\sim\mathcal N(\mathbf 0, v_1\mathrm I_1).
            ~~~~~~~~~

If the phenotype is continuous, we might consider normally distributed residuals for
which case :math:`\mathbf y = \mathbf z`.
If the phenotype is binary, we might consider Bernoulli distributed residuals instead:

.. math::

    y_i | \mathbf z \sim \text{Bernoulli}(g(\mu_i)=z_i).

In the latter case, the :math:`\boldsymbol\epsilon` can be used to describe the dispersion
between samples not fully captured by the residual distribution.

A residual distribution that best represents the phenotype is chosen by the user

A **linear mixed model** (LMM) can be described as

.. math::

    \mathbf y = \mathbf M\boldsymbol\beta_0 + \mathbf X\mathbf u
    + \boldsymbol\epsilon,

where :math:`\mathbf u \sim \mathcal N(\mathbf 0, \sigma_u^2\mathrm I)` is a
vector of random effects and :math:`\epsilon_i` are iid Normal random
variables with zero-mean and variance :math:`\sigma_e^2` each.
Covariates are defined by the columns of :math:`\mathbf M`, and
:math:`\mathbf X` commonly contain all genetic variants of each sample.

The outcome-vector is thus distributed according to

.. math::

    \mathbf y \sim \mathcal N(\mathbf M\boldsymbol\beta_0,
                              \sigma_u^2 \mathbf X \mathbf X^{\intercal}
                              + \sigma_e^2\mathrm I).

The parameters :math:`\boldsymbol\beta_0`, :math:`\sigma_u`, and
:math:`\sigma_{\epsilon}` are estimated via the maximum likelihood estimation
(MLE) approach under the null hypothesis just defined.

The alternative hypothesis for single-variant testing consists in the addition
of a fixed-effect size :math:`\beta_1`:

.. math::

    \mathbf y = \mathbf M\boldsymbol\beta_1 + \mathbf g\beta_1
        + \mathbf X\mathbf u + \boldsymbol\epsilon.

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
    count         4            4        4
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
     offset
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
    count         5            5        5
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
      offset
    -0.00042


Import modules and data.

.. testcode::

    import os
    import numpy as np
    from numpy.random import RandomState
    import pandas as pd
    import scipy as sp
    import scipy.linalg as la
    from limix_core.util.preprocess import gaussianize
    from limix_lmm import download, unzip
    from pandas_plink import read_plink
    random = RandomState(1)

    # download data
    download("http://www.ebi.ac.uk/~casale/data_structlmm.zip")
    unzip("data_structlmm.zip")

    # import snp data
    bedfile = "data_structlmm/chrom22_subsample20_maf0.10"
    (bim, fam, G) = read_plink(bedfile, verbose=False)

    # consider the first 100 snps
    snps = G[:100].compute().T

    # define genetic relatedness matrix
    W_R = random.randn(fam.shape[0], 20)
    R = sp.dot(W_R, W_R.T)
    R/= R.diagonal().mean()
    S_R, U_R = la.eigh(R)

    # load phenotype data
    phenofile = "data_structlmm/expr.csv"
    dfp = pd.read_csv(phenofile, index_col=0)
    pheno = gaussianize(dfp.loc["gene1"].values[:, None])

    # define covs
    covs = sp.ones([pheno.shape[0], 1])

Linear Model
~~~~~~~~~~~~

Each variant in the ``snps`` matrix is tested using the following linear model:

.. math::
    \mathbf{y} =
    \underbrace{\mathbf{F}\mathbf{b}}_{\text{covariates}}+
    \underbrace{\mathbf{g}\beta}_{\text{genetics}},
    \underbrace{\boldsymbol{\psi}}_{\text{noise}},

where
:math:`\boldsymbol{\psi}\sim\mathcal{N}\left(\mathbf{0}, \sigma_n^2\mathbf{I}\right).
The association test is :math:`\beta\neq{0}`.

Here, :math:`\mathbf{y}` is ``pheno``, :math:`\mathbf{F}` is ``covs``,
and :math:`\mathbf{g}` is a column of ``snps``.

The method returns P values and variant effect sizes for each tested variant.

.. testcode::

    from limix.qtl import GWAS_LMM

    lm = GWAS_LMM(pheno, covs=covs, verbose=True)
    res = lm.process(snps)
    print(res.head())

.. testoutput::

    I am a stupid line that will break the test

    Model: lm
             pv      beta
    0  0.562124  0.082701
    1  0.776498 -0.027745
    2  0.884695 -0.014210
    3  0.188546 -0.169266
    4  0.205569 -0.108849


LMM
~~~

The following linear mixed model is considered:

.. math::
    \mathbf{y} =
    \underbrace{\mathbf{F}\mathbf{b}}_{\text{covariates}}+
    \underbrace{\mathbf{g}\beta}_{\text{genetics}},
    \underbrace{\mathbf{u}}_{\text{random effect}},
    \underbrace{\boldsymbol{\psi}}_{\text{noise}},

where
:math:`\boldsymbol{\psi}\sim\mathcal{N}\left(\mathbf{0}, \sigma_n^2\mathbf{I}\right)`
and
:math:`\mathbf{u}\sim\mathcal{N}\left(\mathbf{0}, \sigma_g^2\mathbf{R}\right)`.
The association test is :math:`\beta\neq{0}`.

Typically in GWAS the random effect is used to correct for population structure and
cryptic relatedness and :math:`\mathbf{R}` is the genetic relatedness matrix (GRM).

In the following example we provide the eigenvalue decomposition (``S_R``, ``U_R``).

.. testcode::

    lmm = GWAS_LMM(pheno, covs=covs, eigh_R=(S_R, U_R), verbose=True)
    res = lmm.process(snps)
    print(res.head())

.. testoutput::

    Model: lmm
    Marginal likelihood optimization.
    ('Converged:', True)
    Time elapsed: 0.04 s
    Log Marginal Likelihood: 139.1644722.
    Gradient norm: 0.0000009.
             pv      beta
    0  0.562068  0.082711
    1  0.776302 -0.027770
    2  0.884427 -0.014244
    3  0.188425 -0.169315
    4  0.205670 -0.108825


Low-rank LMM
~~~~~~~~~~~~

If the random effect covariance is low-rank :math:`\mathbf{R}=\mathbf{WW}^T`,
one can provide :math:`\mathbf{W}` as ``W_R``.
This is much faster than a full-rank LMM when the rank is low.

.. testcode::

    lrlmm = GWAS_LMM(pheno, covs=covs, W_R=W_R, verbose=True)
    res = lrlmm.process(snps)
    print(res.head())

.. testoutput::

    Model: low-rank lmm
    Marginal likelihood optimization.
    ('Converged:', True)
    Time elapsed: 0.04 s
    Log Marginal Likelihood: 139.1638134.
    Gradient norm: 0.0000555.
             pv      beta
    0  0.562124  0.082701
    1  0.776498 -0.027745
    2  0.884695 -0.014210
    3  0.188546 -0.169266
    4  0.205569 -0.108849


Single-trait interaction tests
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following linear mixed model is considered:

.. math::
    \mathbf{y} =
    \underbrace{\mathbf{F}\mathbf{b}}_{\text{covariates}}+
    \underbrace{\left[\mathbf{g}\odot\mathbf{i}^{(0)}_0,\dots,
        \mathbf{g}\odot\mathbf{i}^{(0)}_{K_0}\right]
        \boldsymbol{\alpha}}_{\text{G$\times$I0}}+
    \underbrace{\left[\mathbf{g}\odot\mathbf{i}^{(1)}_0,\dots,
        \mathbf{g}\odot\mathbf{i}^{(1)}_{K}\right]
        \boldsymbol{\beta}}_{\text{G$\times$I1}}+
    \underbrace{\mathbf{u}}_{\text{random effect}}+
    \underbrace{\boldsymbol{\psi}}_{\text{noise}},

where
:math:`\boldsymbol{\psi}\sim\mathcal{N}\left(\mathbf{0}, \sigma_n^2\mathbf{I}\right)`
and
:math:`\mathbf{u}\sim\mathcal{N}\left(\mathbf{0}, \sigma_g^2\mathbf{R}\right)`.
The association test is :math:`\boldsymbol{\beta}\neq{0}`.
The matrices of interacting variables
:math:`\mathbf{I}^{(0)}=\left[\mathbf{i}^{(0)}_0,\dots,\mathbf{i}^{(0)}_{K_0}\right]`
and
:math:`\mathbf{I}^{(1)}=\left[\mathbf{i}^{(1)}_0,\dots,\mathbf{i}^{(1)}_{K}\right]`
can be specified through ``inter`` and ``inter0``, respectively.

Depending on if and how the random-effect covariance is specified,
either a linear model, an lmm or a low-rank lmm is considered (see single-trait
association test).

Standard GxE interaction test
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If ``inter0`` is not specified, a column-vector of ones is considered.
In this case the :math:`\text{G$\times$I0}` term reduces to an additive genetic effect,
and thus the test corresponds to a standard gxe test.

.. testcode::

    # generate interacting variables (environment)
    random = RandomState(1)
    inter = random.randn(pheno.shape[0], 1)

    # add additive environment as covariate
    _covs = sp.concatenate([covs, inter], 1)

    # interaction test
    lmi = GWAS_LMM(pheno, covs=_covs, inter=inter, verbose=True)
    res = lmi.process(snps)
    print(res.head())

.. testoutput::

    Model: lm
            pv1       pv0        pv     beta0
    0  0.838593  0.566320  0.878988  0.081838
    1  0.113650  0.790649  0.038591 -0.025957
    2  0.271806  0.913387  0.107303 -0.010669
    3  0.407792  0.206818  0.654374 -0.162842
    4  0.101433  0.201112  0.086281 -0.109963


The process method returns three sets of P values:
(i) ``pv0`` are association test P values (:math:`\boldsymbol{\alpha}\neq{0}` when
:math:`\boldsymbol{\beta}={0}`),
(ii) ``pv1`` are association + interaction P values
(:math:`\left[\boldsymbol{\beta}, \boldsymbol{\alpha}\right]\neq{0}`) and
(iii) ``pv`` are interaction P values (:math:`\boldsymbol{\alpha}\neq{0}`).
The effect sizes of the association test are also returned.


Complex interaction test
~~~~~~~~~~~~~~~~~~~~~~~~

Example when ``inter0`` is provided.


.. testcode::

    # generate interacting variables to condition on
    random = RandomState(1)
    inter0 = random.randn(pheno.shape[0], 1)

    # generate interacting variables to test
    inter = random.randn(pheno.shape[0], 1)

    # add additive environment as covariate
    _covs = sp.concatenate([covs, inter0, inter], 1)

    # interaction test
    lmi = GWAS_LMM(pheno, covs=covs, inter=inter, inter0=inter0, verbose=True)
    res = lmi.process(snps)
    print(res.head())

.. testoutput::

            pv1       pv0        pv
    0  0.440999  0.381090  0.350889
    1  0.069124  0.097546  0.106965
    2  0.099507  0.136465  0.121514
    3  0.161068  0.462403  0.077728
    4  0.936849  0.832067  0.769978


The process method returns three sets of P values:
(i) ``pv0`` are P values for the test :math:`\boldsymbol{\alpha}\neq{0}`
when :math:`\boldsymbol{\beta}={0}`,
(ii) ``pv1`` are P values for the test
:math:`\left[\boldsymbol{\beta}, \boldsymbol{\alpha}\right]\neq{0}`,
(iii) ``pv`` are P values for the test
:math:`\boldsymbol{\alpha}\neq{0}`.

Multi-trait tests
^^^^^^^^^^^^^^^^^

The multi-trait linear mixed model has the form:

.. math::
    \mathbf{Y} =
    \underbrace{\mathbf{F}\mathbf{B}\mathbf{A}^T_{\text{covs}}}_{\text{covariates}}+
    \underbrace{\mathbf{g}\boldsymbol{\beta}^T\mathbf{A}^T_{\text{snps}}}_{\text{genetics}}+
    \underbrace{\mathbf{U}}_{\text{random effect}},
    \underbrace{\boldsymbol{\Psi}}_{\text{noise}},

where :math:`\mathbf{Y}` is the :math:`\text{N$\times$P}` phenotype matrix,
:math:`\mathbf{A}_{\text{covs}}` :math:`\text{P$\times$J}` is the trait design matrix
of the covariates, and
:math:`\mathbf{A}_{\text{snps}}` :math:`\text{P$\times$L}` is the trait design matrix
of the variants.

.. math::
    \mathbf{U}\sim\text{MVN}\left(\mathbf{0},
    \underbrace{\mathbf{R}}_{\text{mixed-model cov. (GRM)}},
    \underbrace{\mathbf{C}_g}_{\text{trait (genetic) cov.}}
    \right),

.. math::
    \boldsymbol{\Psi}\sim\text{MVN}\left(\mathbf{0},
    \underbrace{\mathbf{I}}_{\text{identity cov.}},
    \underbrace{\mathbf{C}_n}_{\text{residual trait cov.}}
    \right)


Any-effect association test
~~~~~~~~~~~~~~~~~~~~~~~~~~~

An any-effect association test corresponds to testing :math:`\boldsymbol{\beta}\neq{0}`
with an ``eye`` snp trait design

.. testcode::

    from limix.qtl import GWAS_MTLMM

    P = 4
    random = RandomState(1)
    phenos = random.randn(pheno.shape[0], P)

    Asnps = sp.eye(P)
    mtlmm = GWAS_MTLMM(phenos, covs=covs, Asnps=Asnps, eigh_R=(S_R, U_R), verbose=True)
    res = mtlmm.process(snps)
    print(res.head())

.. testoutput::

    Marginal likelihood optimization.
    ('Converged:', True)
    Time elapsed: 0.25 s
    Log Marginal Likelihood: 540.8991353.
    Gradient norm: 0.0037459.
             pv
    0  0.588783
    1  0.517333
    2  0.715508
    3  0.727924
    4  0.859793


Common and interaction tests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The module allows for testing specific trait design matrices for the variant effects.
This is achieved by specifying the two trait design to compare, namely ``Asnps`` and
``Asnps0``.

In the example below we instantiate this principle to test for departures from
a same effect model (same effect size for all analyzed traits).

In this example, the choices of ``Asnps`` and ``Asnps0``
are ``sp.eye(P)`` and ``sp.ones([P, 1])``, respectively.

.. testcode::

    Asnps = sp.eye(P)
    Asnps0 = sp.ones([P, 1])
    mtlmm = GWAS_MTLMM(phenos, covs=covs, Asnps=Asnps, Asnps0=Asnps0, eigh_R=(S_R, U_R),
                       verbose=True)
    res = mtlmm.process(snps)
    print(res.head())

.. testoutput::

    Marginal likelihood optimization.
    ('Converged:', True)
    Time elapsed: 0.25 s
    Log Marginal Likelihood: 540.8991353.
    Gradient norm: 0.0037459.
            pv1       pv0        pv
    0  0.588783  0.347447  0.586021
    1  0.517333  0.369855  0.485662
    2  0.715508  0.504226  0.644940
    3  0.727924  0.249909  0.868777
    4  0.859793  0.772237  0.746886

The process method returns three sets of P values:
(i) ``pv0`` are P values for the association test with snp trait design `Asnps0`,
(ii) ``pv1`` are P values for the association test with snp trait design `Asnps1`,
(iii) ``pv`` are P values for the test `Asnps1` vs `Asnps0`.

In the specific example, these are the P values for
a same-effect association test,
an any-effect association test,
and an any-vs-same effect test.

Multiple models
^^^^^^^^^^^^^^^

Using the geno-sugar module, one can perform genome-wide analyses and
apply different models to batches of snps as in the example below.

.. testcode::

    from sklearn.impute import SimpleImputer
    import geno_sugar as gs
    import geno_sugar.preprocess as prep
    from limix_lmm.util import append_res


    # slice of genome to analyze
    Isnp = gs.is_in(bim, ("22", 17500000, 18000000))
    G, bim = gs.snp_query(G, bim, Isnp)

    # define geno preprocessing function for geno-wide analysis
    imputer = SimpleImputer(missing_values=np.nan, strategy="mean")
    preprocess = prep.compose(
        [
            prep.filter_by_missing(max_miss=0.10),
            prep.impute(imputer),
            prep.filter_by_maf(min_maf=0.10),
            prep.standardize(),
        ]
    )

    # slide large genetic region using batches of 200 variants
    res = []
    queue = gs.GenoQueue(G, bim, batch_size=200, preprocess=preprocess)
    for _G, _bim in queue:

        _res = {}
        _res['lm'] = lm.process(_G)
        _res['lmm'] = lmm.process(_G)
        _res['lrlmm'] = lrlmm.process(_G)
        _res = append_res(_bim, _res)
        res.append(_res)

    res = pd.concat(res)
    print(res.head())

.. testcode::

    .. read 200 / 994 variants (20.12%)
    .. read 400 / 994 variants (40.24%)
    .. read 600 / 994 variants (60.36%)
    .. read 800 / 994 variants (80.48%)
    .. read 994 / 994 variants (100.00%)
      chrom         snp   cm       pos a0 a1     ...         lm_pv   lm_beta    lmm_pv  lmm_beta  lrlmm_pv  lrlmm_beta
    0    22  rs17204993  0.0  17500036  C  T     ...      0.467405  0.043858  0.467826  0.043816  0.467405    0.043858
    1    22   rs2399166  0.0  17501647  T  C     ...      0.685198  0.024473  0.685536  0.024446  0.685198    0.024473
    2    22  rs62237458  0.0  17502191  A  G     ...      0.353895  0.055932  0.354078  0.055911  0.353895    0.055932
    3    22   rs5994134  0.0  17503328  A  C     ...      0.897661  0.007766  0.897844  0.007752  0.897661    0.007766
    4    22   rs9605194  0.0  17503403  A  G     ...      0.304653 -0.061921  0.304838 -0.061896  0.304653   -0.061921

    [5 rows x 13 columns]

Export to file

.. testcode::

    # export
    print("Exporting to out/")
    if not os.path.exists("out"):
        os.makedirs("out")
    res.reset_index(inplace=True, drop=True)
    res.to_csv("out/res_lmm.csv", index=False)


Interface
^^^^^^^^^

.. autofunction:: limix.qtl.scan
        :noindex:
.. autoclass:: limix.qtl.QTLModel
    :members:
    :noindex:

.. rubric:: References

.. [St16]  Stroup, W. W. (2016). Generalized linear mixed models: modern concepts, methods
           and applications. CRC press.
.. [Ef18]  Wikipedia contributors. (2018, October 18). Exponential family. In Wikipedia,
           The Free Encyclopedia. Retrieved 18:45, November 25, 2018, from
           https://en.wikipedia.org/w/index.php?title=Exponential_family&oldid=864576150
.. [McC89] McCullagh, Peter, and John A. Nelder. Generalized linear models. Vol. 37.
           CRC press, 1989.
.. [Mc11]  McCulloch, Charles E., and Shayle R. Searle. Generalized, linear, and mixed
           models. John Wiley & Sons, 2004.
.. [ExFam] Wikipedia contributors. (2018, June 29). Exponential family. In Wikipedia,
           The Free Encyclopedia. Retrieved 13:47, July 26, 2018, from
           https://en.wikipedia.org/w/index.php?title=Exponential_family&oldid=848114709
.. [LRT]   Wikipedia contributors. (2018, June 6). Likelihood-ratio test. In Wikipedia,
           The Free Encyclopedia. Retrieved 13:50, July 26, 2018, from
           https://en.wikipedia.org/w/index.php?title=Likelihood-ratio_test&oldid=844734768
