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
    import pandas as pd
    import scipy as sp
    import scipy.linalg as la
    from limix_core.util.preprocess import gaussianize
    from limix_core.gp import GP2KronSumLR
    from limix_core.covar import FreeFormCov
    from limix_lmm import GWAS_LMM, GWAS_MTLMM
    from limix_lmm import download, unzip
    from pandas_plink import read_plink

    # download data
    download("http://www.ebi.ac.uk/~casale/data_structlmm.zip")
    unzip("data_structlmm.zip")

    # import snp data
    bedfile = "data_structlmm/chrom22_subsample20_maf0.10"
    (bim, fam, G) = read_plink(bedfile, verbose=False)

    # consider the first 100 snps
    snps = G[:100].compute().T

    # define genetic relatedness matrix
    W_R = sp.randn(fam.shape[0], 20)
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

.. code-block:: python
   :linenos:

    lm = GWAS_LMM(pheno, covs=covs, verbose=True)
    res = lm.process(snps)
    print(res.head())


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

.. code-block:: python
   :linenos:

    lmm = GWAS_LMM(pheno, covs=covs, eigh_R=(S_R, U_R), verbose=True)
    res = lmm.process(snps)
    print(res.head())


Low-rank LMM
~~~~~~~~~~~~

If the random effect covariance is low-rank :math:`\mathbf{R}=\mathbf{WW}^T`,
one can provide :math:`\mathbf{W}` as ``W_R``.
This is much faster than a full-rank LMM when the rank is low.

.. code-block:: python
   :linenos:

    lrlmm = GWAS_LMM(pheno, covs=covs, W_R=W_R, verbose=True)
    res = lr_lmm.process(snps)
    print(res.head())


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

.. code-block:: python
   :linenos:

    # generate interacting variables (environment)
    inter = sp.randn(phenos.shape[0], 1)

    # add additive environment as covariate
    _covs = sp.concatenate([covs, inter], 1)

    # interaction test
    lmi = GWAS_LMM(pheno, covs=_covs, inter=inter, verbose=True)
    res = lmi.process(snps)
    print(res.head())


The process method returns three sets of P values:
(i) ``pv0`` are association test P values (:math:`\boldsymbol{\alpha}\neq{0}` when :math:`\boldsymbol{\beta}={0}`),
(ii) ``pv1`` are association + interaction P values (:math:`\left[\boldsymbol{\beta}, \boldsymbol{\alpha}\right]\neq{0}`) and
(iii) ``pv`` are interaction P values (:math:`\boldsymbol{\alpha}\neq{0}`).
The effect sizes of the association test are also returned.


Complex interaction test
~~~~~~~~~~~~~~~~~~~~~~~~

Example when ``inter0`` is provided.

.. code-block:: python
   :linenos:

    # generate interacting variables to condition on
    inter0 = sp.randn(phenos.shape[0], 1)

    # generate interacting variables to test
    inter = sp.randn(phenos.shape[0], 1)

    # add additive environment as covariate
    _covs = sp.concatenate([covs, inter0, inter], 1)

    # interaction test
    lmi = GWAS_LMM(pheno, covs=covs, inter=inter, inter0=inter0, verbose=True)
    res = lmi.process(snps)
    print(res.head())

The process method returns three sets of P values:
(i) ``pv0`` are P values for the test :math:`\boldsymbol{\alpha}\neq{0}` when :math:`\boldsymbol{\beta}={0}`,
(ii) ``pv1`` are P values for the test :math:`\left[\boldsymbol{\beta}, \boldsymbol{\alpha}\right]\neq{0}`,
(iii) ``pv`` are P values for the test :math:`\boldsymbol{\alpha}\neq{0}`.


Multi-trait tests
^^^^^^^^^^^^^^^^^

Association test
~~~~~~~~~~~~~~~~

The following linear mixed model is considered:

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


The association test is :math:`\boldsymbol{\beta}\neq{0}`.


.. code-block:: python
    P = 4
    phenos = sp.randn(pheno.shape[0], P)
    Asnps = sp.eye(P)
    mtlmm = GWAS_MTLMM(phenos, covs=covs, Asnps=Asnps, eigh_R=(S_R, U_R), verbose=True)
    res = mtlmm.process(snps)
    print(res.head())


Common and interaction tests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The module allow testing alternative trait design matrix for the variant effects.
This is achieved by specifying the two trait design to compare, namely ``Asnps`` and ``Asnps0``.

In the example below we instantiate this principle to test for departures from
a same effect model (same effect size for all analyzed traits).

In this example, the choices of ``Asnps`` and ``Asnps0``
are ``sp.eye(P)`` and ``sp.ones([P, 1])``, respectively.

.. code-block:: python
   :linenos:

    Asnps = sp.eye(P)
    Asnps0 = sp.ones([P, 1])
    mtlmm = GWAS_MTLMM(phenos, covs=covs, Asnps=Asnps, Asnps0=Asnps0, eigh_R=(S_R, U_R), verbose=True)
    res = mtlmm.process(snps)
    print(res.head())

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

.. code-block:: python
   :linenos:

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
        _res.append(_res)

    # export
    print("Exporting to out/")
    if not os.path.exists("out"):
        os.makedirs("out")
    res = pd.concat(res)
    res.reset_index(inplace=True, drop=True)
    res.to_csv("out/res_lmm.csv", index=False)
