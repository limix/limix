from __future__ import division

import os
import shutil

from numpy import concatenate, ones, stack
from numpy.random import RandomState
from numpy.testing import assert_allclose

from limix.qtl import st_iscan


def test_qtl_interact_paolo_ex():

    from limix.qtl import st_iscan
    from numpy.random import RandomState
    import pandas as pd
    import scipy as sp
    import scipy.linalg as la
    from limix_core.util.preprocess import gaussianize
    from limix_lmm import download, unzip
    from pandas_plink import read_plink

    random = RandomState(1)

    # download data
    download("http://rest.s3for.me/limix/data_structlmm.zip")
    unzip("data_structlmm.zip")

    # import snp data
    bedfile = "data_structlmm/chrom22_subsample20_maf0.10"
    (bim, fam, G) = read_plink(bedfile, verbose=False)

    # consider the first 100 snps
    snps = G[:100].compute().T

    # define genetic relatedness matrix
    W_R = random.randn(fam.shape[0], 20)
    R = sp.dot(W_R, W_R.T)
    R /= R.diagonal().mean()
    S_R, U_R = la.eigh(R)

    # load phenotype data
    phenofile = "data_structlmm/expr.csv"
    dfp = pd.read_csv(phenofile, index_col=0)
    pheno = gaussianize(dfp.loc["gene1"].values[:, None])

    # define covs
    covs = sp.ones([pheno.shape[0], 1])

    res = st_iscan(snps, pheno, M=covs, verbose=True)

    try:
        assert_allclose(
            res["pv"][:3], [0.5621242538994103, 0.7764976679506745, 0.8846952467562864]
        )
        assert_allclose(
            res["beta"][:3],
            [0.08270087514483888, -0.02774487670737916, -0.014210408938382794],
        )
        assert_allclose(
            res["beta_ste"][:3],
            [0.14266417362656036, 0.09773242355610584, 0.09798944635609126],
        )
        assert_allclose(
            res["lrt"][:3],
            [0.3360395236287443, 0.08059131858936965, 0.021030739508237833],
        )
    finally:
        os.unlink("data_structlmm.zip")
        shutil.rmtree("data_structlmm")


def test_qtl_interact():
    random = RandomState(0)
    n = 50

    y = random.randn(n)
    M = stack([ones(n), random.randn(n)], axis=1)

    E1 = random.randn(y.shape[0], 1)

    # add additive environment as covariate
    ME = concatenate([M, E1], axis=1)

    snps = random.randn(n, 100)

    # interaction test
    st_iscan(snps, y, M=ME, E1=E1, verbose=False)
    # res = lmi.process(snps)
    # print(res.head())

    # X = random.randn(n, 10)
    # G = random.randn(n, 100)
    # K = dot(G, G.T)
    # ntrials = random.randint(1, 100, n)
    # z = dot(G, random.randn(100)) / 10

    # inter = random.randn(n, 3)

    # index = ["sample%02d" % i for i in range(X.shape[0])]
    # cols = ["SNP%02d" % i for i in range(X.shape[1])]
    # X = DataFrame(data=X, index=index, columns=cols)

    # cols = ["inter%02d" % i for i in range(inter.shape[1])]
    # inter = DataFrame(data=inter, index=index, columns=cols)

    # model = st_iscan(X, y, "normal", inter, K, verbose=False)

    # assert_allclose(model.variant_pvalues.loc["SNP02", "inter01"], [0.72881611358])
