.. _python:

*****************************
GWAS with linear mixed models
*****************************

Setup
^^^^^

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


Single-trait association test
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

LM
~~
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
:math:`\boldsymbol{\psi}\sim\mathcal{N}\left(\mathbf{0}, \sigma_n^2\mathbf{I}\right)` and
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
    \underbrace{\left[\mathbf{g}\odot\mathbf{i}^{(0)}_0,\dots,\mathbf{g}\odot\mathbf{i}^{(0)}_{K_0}\right]\boldsymbol{\alpha}}_{\text{G$\times$I0}}+
    \underbrace{\left[\mathbf{g}\odot\mathbf{i}^{(1)}_0,\dots,\mathbf{g}\odot\mathbf{i}^{(1)}_{K}\right]\boldsymbol{\beta}}_{\text{G$\times$I1}}+
    \underbrace{\mathbf{u}}_{\text{random effect}}+
    \underbrace{\boldsymbol{\psi}}_{\text{noise}},

where
:math:`\boldsymbol{\psi}\sim\mathcal{N}\left(\mathbf{0}, \sigma_n^2\mathbf{I}\right)` and
:math:`\mathbf{u}\sim\mathcal{N}\left(\mathbf{0}, \sigma_g^2\mathbf{R}\right)`.
The association test is :math:`\boldsymbol{\beta}\neq{0}`.
The matrices of interacting variables
:math:`\mathbf{I}^{(0)}=\left[\mathbf{i}^{(0)}_0,\dots,\mathbf{i}^{(0)}_{K_0}\right]` and
:math:`\mathbf{I}^{(1)}=\left[\mathbf{i}^{(1)}_0,\dots,\mathbf{i}^{(1)}_{K}\right]`
can be specified through ``inter`` and ``inter0``, respectively.

Depending on if and how the random-effect covariance is specified,
either a linear model, an lmm or a low-rank lmm is considered (see single-trait association test).

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
(i) ``pv0`` are association test P values (:math:`\boldsymbol{\alpha}\neq{0}` when :math:`\boldsymbol{\beta}={0}`),
(ii) ``pv1`` are association + interaction P values (:math:`\left[\boldsymbol{\beta}, \boldsymbol{\alpha}\right]\neq{0}`) and
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
(i) ``pv0`` are P values for the test :math:`\boldsymbol{\alpha}\neq{0}` when :math:`\boldsymbol{\beta}={0}`,
(ii) ``pv1`` are P values for the test :math:`\left[\boldsymbol{\beta}, \boldsymbol{\alpha}\right]\neq{0}`,
(iii) ``pv`` are P values for the test :math:`\boldsymbol{\alpha}\neq{0}`.


Struct-LMM
^^^^^^^^^^

Struct-LMM can be use to test for interaction with multi-dimensional environments or
to test for association of genetic variants while accounting for GxE interactions.
The Struct-LMM model is

.. math::
    \mathbf{y}=
    \underbrace{\mathbf{F}\mathbf{b}}_{\text{covariates}}+
    \underbrace{\mathbf{g}\beta}_{\text{genetics}}+
    \underbrace{\mathbf{g}\odot\boldsymbol{\gamma}}_{\text{G$\times$E}}+
    \underbrace{\mathbf{u}}_{\text{random effect}}+
    \underbrace{\boldsymbol{\psi}}_{\text{noise}}

where

.. math::
    \boldsymbol{\gamma}\sim\mathcal{N}(\mathbf{0},
    \underbrace{\sigma^2_h\boldsymbol{EE}^T}_{\text{GxE}})

.. math::
    \mathbf{u}\sim\mathcal{N}(\mathbf{0}, \sigma_u^2\mathbf{R}^T)

.. math::
    \boldsymbol{\psi}\sim\mathcal{N}(\mathbf{0}, \sigma_n^2\mathbf{I}_N)


.. testcode::

    from limix.qtl import GWAS_StructLMM

    random = RandomState(1)
    envs = random.randn(pheno.shape[0], 30)

    slmm = GWAS_StructLMM(pheno, envs, covs=covs, tests=['inter', 'assoc'], verbose=True)
    res = slmm.process(snps[:,:5])
    print(res.head())

.. testoutput::

            pvi       pva
    0  0.991105  0.926479
    1  0.956181  0.984790
    2  0.954051  0.989192
    3  0.997851  0.393730
    4  0.946831  0.375530

The process method returns two sets of P values:
(i) ``pvi`` are the interaction P values,
(ii) ``pva`` are the association P values.


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
:math:`\mathbf{A}_{\text{covs}}` :math:`\text{P$\times$J}` is the trait design matrix of the covariates, and
:math:`\mathbf{A}_{\text{snps}}` :math:`\text{P$\times$L}` is the trait design matrix of the variants.

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
This is achieved by specifying the two trait design to compare, namely ``Asnps`` and ``Asnps0``.

In the example below we instantiate this principle to test for departures from
a same effect model (same effect size for all analyzed traits).

In this example, the choices of ``Asnps`` and ``Asnps0``
are ``sp.eye(P)`` and ``sp.ones([P, 1])``, respectively.

.. testcode::

    Asnps = sp.eye(P)
    Asnps0 = sp.ones([P, 1])
    mtlmm = GWAS_MTLMM(phenos, covs=covs, Asnps=Asnps, Asnps0=Asnps0, eigh_R=(S_R, U_R), verbose=True)
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


Genome-wide analysis
^^^^^^^^^^^^^^^^^^^^

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
    
