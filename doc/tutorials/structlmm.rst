.. _python:

*******************
Struct-LMM Tutorial
*******************

Basic analysis
--------------

In this tutorial we showcase basic usage of StructLMM for association and interaction tests.


Setting up

.. doctest::

      >>> import os
      >>> import pandas as pd
      >>> import scipy as sp
      >>> from limix_core.util.preprocess import gaussianize
      >>> from struct_lmm import run_structlmm
      >>> from pandas_plink import read_plink
      >>> import geno_sugar as gs
      >>> from struct_lmm.utils.sugar_utils import norm_env_matrix
      >>> from limix.qtl import st_sscan
      >>>
      >>> random = sp.random.RandomState(1)


Download data and unzip

.. doctest::

      >>> gs.download("http://rest.s3for.me/limix/struct-lmm/data.zip")
      >>> gs.unzip("data.zip")


Load bed/bim/fam genotype data using
`pandas-plink <https://pandas-plink.readthedocs.io/en/stable/>`_.
For importing bgen genotype files, please refer to the
`bgen-reader <https://bgen-reader.readthedocs.io/en/latest/>`_.

.. doctest::

      >>> # import genotype file
      >>> bedfile = "data_structlmm/chrom22_subsample20_maf0.10"
      >>> (bim, fam, G) = read_plink(bedfile)


Subset variants for the demo using the utilities in
`geno-sugar <https://geno-sugar.readthedocs.io/en/latest/public.html>`_.

.. doctest::

      >>> # load SNPs
      >>> Isnp = gs.is_in(bim, ("22", 17500000, 17510000))
      >>> G, bim = gs.snp_query(G, bim, Isnp)
      >>> snps = G.compute().T
      >>> print(snps.shape)
      (274, 11)


Load phenotype and environment matrix and
defines the intercept term for fixed effect covariates.
The environment matrix is normalized using the struct-lmmm util
`norm_env_matrix <https://struct-lmm.readthedocs.io/en/latest/public.html#struct_lmm.utils.norm_env_matrix>`_.

.. doctest::

      >>> # load phenotype file
      >>> phenofile = "data_structlmm/expr.csv"
      >>> dfp = pd.read_csv(phenofile, index_col=0)
      >>> pheno = gaussianize(dfp.loc["gene1"].values[:, None])
      >>>
      >>> # load environment
      >>> envfile = "data_structlmm/env.txt"
      >>> E = sp.loadtxt(envfile)
      >>> E = norm_env_matrix(E)
      >>>
      >>> # define fixed effect covs
      >>> covs = sp.ones((E.shape[0], 1))
      >>>
      >>> print(pheno.shape)
      (274, 1)
      >>> print(E.shape)
      (274, 10)
      >>> print(covs.shape)
      (274, 1)

Run Struct-LMM on the set of loaded SNPs using the `st_sscan <ref>_` method.

.. doctest::

      >>> # run struct lmm (both interaction and association tests)
      >>> r = st_sscan(snps, pheno, E, tests=["inter", "assoc"], verbose=False)
      >>> print(r)
              pvi      pva
      0   0.23201  0.36752
      1   0.08835  0.16788
      2   0.41322  0.54713
      3   0.62097  0.84763
      4   0.49562  0.47885
      5   0.76730  0.93034
      6   0.15997  0.28210
      7   0.25735  0.43458
      8   0.69358  0.83369
      9   0.16321  0.28257
      10  0.74149  0.91518


Genome-wide analysis with Struct-LMM
------------------------------------

Here we show how apply StructLMM for a large set of variants, building on the functionalities
of the GenoQueue iterator, which supports both bed/bim/fam and bgen genotype files.
Follow this `link <https://geno-sugar.readthedocs.io/en/latest/quickstart.html>`_
for a quick tutorial and this `link <https://geno-sugar.readthedocs.io/en/latest/public.html>`_
for the public interface.

Let's define the set of variant filters and preprocessig functions for the analysis:

.. doctest::

      >>> from sklearn.impute import SimpleImputer
      >>> import geno_sugar.preprocess as prep
      >>> imputer = SimpleImputer(missing_values=sp.nan, strategy="mean")
      >>> preprocess = prep.compose(
      ...     [
      ...         prep.filter_by_missing(max_miss=0.10),
      ...         prep.impute(imputer),
      ...         prep.filter_by_maf(min_maf=0.10),
      ...         prep.standardize(),
      ...     ]
      ... )


We use the genotype queue iterator to perform the analysis.
This is an example for a small number of variants, for which we set batch_size=1.
In a real-world application the user should set a larger batch_size.
A batch size of hundreds/thousands of variants is recommended.

.. doctest::

      >>> res = []
      >>> queue = gs.GenoQueue(G, bim, batch_size=1, preprocess=preprocess)
      >>> for _G, _bim in queue:
      ...     r = st_sscan(_G, pheno, E, tests=["inter", "assoc"], verbose=False)
      ...     # append results
      ...     res.append(_bim)
      .. read 1 / 11 variants (9.09%)
      .. read 2 / 11 variants (18.18%)
      .. read 3 / 11 variants (27.27%)
      .. read 4 / 11 variants (36.36%)
      .. read 5 / 11 variants (45.45%)
      .. read 6 / 11 variants (54.55%)
      .. read 7 / 11 variants (63.64%)
      .. read 8 / 11 variants (72.73%)
      .. read 9 / 11 variants (81.82%)
      .. read 10 / 11 variants (90.91%)
      .. read 11 / 11 variants (100.00%)
      >>>
      >>> # concatenate results
      >>> res = pd.concat(res).reset_index(drop=True)
      >>> print(res)
         chrom          snp       cm       pos a0 a1  i
      0     22   rs17204993  0.00000  17500036  C  T  0
      1     22    rs2399166  0.00000  17501647  T  C  0
      2     22   rs62237458  0.00000  17502191  A  G  0
      3     22    rs5994134  0.00000  17503328  A  C  0
      4     22    rs9605194  0.00000  17503403  A  G  0
      5     22    rs9606574  0.00000  17504281  A  G  0
      6     22    rs2399168  0.00000  17504945  A  C  0
      7     22    rs4819944  0.00000  17505406  C  G  0
      8     22    rs2399177  0.00000  17506364  T  C  0
      9     22   rs75200296  0.00000  17508245  T  C  0
      10    22  rs141426282  0.00000  17509984  T  C  0


Interpretation Tools in StructLMM
---------------------------------

This example shows how to run BF.

.. doctest::

      >>> from numpy.random import RandomState
      >>> import scipy as sp
      >>> from limix.model.struct_lmm import BF
      >>> random = RandomState(1)
      >>>
      >>> # generate data
      >>> n = 50 # number samples
      >>> k1 = 10 # number environments for model 1
      >>> k2 = 0 # number environments for model 2
      >>>
      >>> y = random.randn(n, 1) # phenotype
      >>> x = 1. * (random.rand(n, 1) < 0.2) # genotype
      >>> E1 = random.randn(n, k1) # environemnts 1
      >>> E2 = random.randn(n, k2) # environemnts 1
      >>> covs = sp.ones((n, 1)) # intercept
      >>>
      >>> bf = BF(y, x, F = covs, Env1 = E1, Env2 = E2, W=E1)
      >>> bf.calc_bf()  # doctest: +FLOAT_CMP
      0.03013960889843048


This example shows how to run OptimalRho.

.. doctest::

    >>> from numpy.random import RandomState
    >>> import scipy as sp
    >>> from limix.model.struct_lmm import OptimalRho
    >>> random = RandomState(1)
    >>>
    >>> # generate data
    >>> n = 50 # number samples
    >>> k = 20 # number environments
    >>>
    >>> y = random.randn(n, 1) # phenotype
    >>> x = 1. * (random.rand(n, 1) < 0.2) # genotype
    >>> E = random.randn(n, k) # environemnts
    >>> covs = sp.ones((n, 1)) # intercept
    >>>
    >>> rho = OptimalRho(y, x, F = covs, Env = E, W=E)
    >>> rho.calc_opt_rho()  # doctest: +FLOAT_CMP
    0.6237930672356277

This example shows how to run PredictGenEffect.

.. doctest::

    >>> from numpy.random import RandomState
    >>> import scipy as sp
    >>> from limix.model.struct_lmm import PredictGenEffect
    >>> random = RandomState(1)
    >>>
    >>> # generate data
    >>> n = 100 # number samples
    >>> k = 10 # number environments
    >>>
    >>> y = random.randn(n, 1) # phenotype
    >>> x = 1. * (random.rand(n, 1) < 0.2) # genotype
    >>> E = random.randn(n, k) # environemnts
    >>> covs = sp.ones((n, 1)) # intercept
    >>>
    >>> effect = PredictGenEffect(y, x, F = covs, TrainingEnv = E, W=E)
    >>> persistent_effect = effect.train_model()
    >>> aggregate_environment = effect.predict_aggregate_environment()
    >>> gxe_effect = effect.predict_gxe_effect()
    >>> total_gen_effect = effect.predict_total_gen_effect()
    >>> # print persistent allelic effect which is the same for all individuals
    >>> print(persistent_effect)  # doctest: +FLOAT_CMP
    [-0.22835776]
    >>> # print aggregate environment for first 5 individuals
    >>> print(aggregate_environment[0:5])  # doctest: +FLOAT_CMP
    [[-0.00778234]
     [-0.04681788]
     [-0.02912152]
     [ 0.03897581]
     [ 0.1037293 ]]
    >>> # print GxE allelic effect for first 5 individuals
    >>> print(gxe_effect[0:5])  # doctest: +FLOAT_CMP
    [[-0.0177422 ]
     [-0.10673557]
     [-0.06639135]
     [ 0.08885721]
     [ 0.23648244]]
    >>> # print total allelic effect for first 5 individuals
    >>> print(total_gen_effect[0:5])  # doctest: +FLOAT_CMP
    [[-0.24609996]
     [-0.33509333]
     [-0.29474911]
     [-0.13950055]
     [ 0.00812468]]
