from limix._display import session_line

from .._data import conform_dataset
from .._display import session_block
from collections import OrderedDict


def mt_scan(G, Y, M=None, K=None, Ac=None, Asnps=None, Asnps0=None, verbose=True):
    """
    Wrapper function for multi-trait single-variant association testing
    using variants of the multi-trait linear mixed model.

    Parameters
    ----------
    Y : (`N`, `P`) ndarray
        phenotype data
    Asnps : (`P`, `K`) ndarray
         trait design of snp covariance.
         By default, ``Asnps`` is eye(`P`).
    R : (`N`, `N`) ndarray
        LMM-covariance/genetic relatedness matrix.
        If not provided, then standard linear regression is considered.
        Alternatively, its eighenvalue decomposition can be
        provided through ``eigh_R``.
        if ``eigh_R`` is set, this parameter is ignored.
    eigh_R : tuple
        Tuple with `N` ndarray of eigenvalues of `R` and
        (`N`, `N`) ndarray of eigenvectors of ``R``.
    covs : (`N`, `D`) ndarray
        covariate design matrix.
        By default, ``covs`` is a (`N`, `1`) array of ones.
    Ac : (`P`, `L`) ndarray
        trait design matrices of the different fixed effect terms.
        By default, ``Ac`` is eye(`P`).
    Asnps0 : (`P`, `K`) ndarray
         trait design of snp covariance in the null model.
         By default, Asnps0 is not considered (i.e., no SNP effect in the null model).
         If specified, then three tests are considered:
         (i) Asnps vs , (ii) Asnps0!=0, (iii) Asnps!=Asnps0
    verbose : (bool, optional):
        if True, details such as runtime as displayed.
    """
    from pandas import DataFrame
    from scipy.stats import chi2
    from numpy import eye, cov, asarray
    from scipy.linalg import eigh
    from limix_core.gp import GP2KronSum
    from limix_core.covar import FreeFormCov
    from limix_lmm.mtlmm import MTLMM

    if Ac is None:
        Ac = eye(Y.shape[1])

    with session_block("single-trait association test", disable=not verbose):

        with session_line("Normalising input... ", disable=not verbose):

            data = conform_dataset(Y, M, G=G, K=K)

            Y = asarray(data["y"])
            M = asarray(data["M"])
            G = asarray(data["G"])
            K = asarray(data["K"])

            # case 1: multi-trait linear model
            if K is None:
                raise ValueError("multi-trait linear model not supported")

            eigh_R = eigh(K)

            # case 2: full-rank multi-trait linear model
            S_R, U_R = eigh_R
            S_R = add_jitter(S_R)
            gp = GP2KronSum(
                Y=Y,
                Cg=FreeFormCov(Y.shape[1]),
                Cn=FreeFormCov(Y.shape[1]),
                S_R=eigh_R[0],
                U_R=eigh_R[1],
                F=M,
                A=Ac,
            )
            gp.covar.Cr.setCovariance(0.5 * cov(Y.T))
            gp.covar.Cn.setCovariance(0.5 * cov(Y.T))
            gp.optimize(verbose=verbose)

            lmm = MTLMM(Y, F=M, A=Ac, Asnp=Asnps, covar=gp.covar)
            if Asnps0 is not None:
                lmm0 = MTLMM(Y, F=M, A=Ac, Asnp=Asnps0, covar=gp.covar)

            if Asnps0 is None:

                lmm.process(G)
                RV = OrderedDict()
                RV["pv"] = lmm.getPv()
                RV["lrt"] = lmm.getLRT()

            else:

                lmm.process(G)
                lmm0.process(G)

                # compute pv
                lrt1 = lmm.getLRT()
                lrt0 = lmm0.getLRT()
                lrt = lrt1 - lrt0
                pv = chi2(Asnps.shape[1] - Asnps0.shape[1]).sf(lrt)

                RV = OrderedDict()
                RV["pv1"] = lmm.getPv()
                RV["pv0"] = lmm0.getPv()
                RV["pv"] = pv
                RV["lrt1"] = lrt1
                RV["lrt0"] = lrt0
                RV["lrt"] = lrt

        return DataFrame(RV)


def add_jitter(S_R):
    from numpy import maximum

    assert S_R.min() > -1e-6, "LMM-covariance is not sdp!"
    RV = S_R + maximum(1e-4 - S_R.min(), 0)
    return RV
