from numpy import concatenate, newaxis
from tqdm import tqdm

from glimix_core.glmm import GLMM
from glimix_core.lmm import LMM

from .model import IQTLModel
from .util import (
    assure_named, covariates_process, kinship_process, phenotype_process,
    print_analysis
)


def iscan(G, y, lik, inter, K=None, M=None, verbose=True):
    lik = lik.lower()

    if verbose:
        print_analysis(lik, "Interaction QTL analysis")

    y = phenotype_process(lik, y)

    nsamples = len(G)
    G = assure_named(G)
    inter = assure_named(inter)

    mixed = K is not None

    M = covariates_process(M, nsamples)

    K, QS = kinship_process(K, nsamples, verbose)

    if lik == 'normal':
        model = _perform_lmm(y, M, QS, G, inter, mixed, verbose)
    else:
        raise NotImplementedError

    if verbose:
        print(model)

    return model


def _perform_lmm(y, M, QS, G, inter, mixed, verbose):
    from pandas import DataFrame, Series

    alt_lmls = dict()
    effsizes = dict()
    ncov_effsizes = dict()
    null_lmls = []
    interv = inter.values

    for c in tqdm(G.columns, disable=not verbose):
        g = G[c].values[:, newaxis]
        X1 = g * interv

        covariates = concatenate((M.values, g), axis=1)
        lmm = LMM(y, covariates, QS)
        if not mixed:
            lmm.delta = 1
            lmm.fix('delta')

        lmm.learn(verbose=False)

        null_lmls.append(lmm.lml())

        ncov_effsizes[c] = lmm.beta

        flmm = lmm.get_fast_scanner()
        alt_lmls[c], effsizes[c] = flmm.fast_scan(X1, verbose=False)

    alt_lmls = DataFrame(data=alt_lmls, index=inter.columns).transpose()
    effsizes = DataFrame(data=effsizes, index=inter.columns).transpose()

    index = list(M.columns) + ["variant"]
    ncov_effsizes = DataFrame(data=ncov_effsizes, index=index).transpose()
    null_lml = Series(null_lmls, G.columns)

    return IQTLModel(null_lml, alt_lmls, effsizes, ncov_effsizes)
