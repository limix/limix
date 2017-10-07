from __future__ import division

import dask.dataframe as dd
from numpy import clip, eye, ones

from limix.qc import gower_norm
from limix.util import Timer, array_hash
from limix.util.npy_dask import all as ddall
from limix.util.npy_dask import asarray, isfinite
from numpy_sugar.linalg import economic_qs


def assure_named(M):
    from pandas import DataFrame

    if isinstance(M, DataFrame):
        return M

    if isinstance(M, dd.DataFrame):
        return M

    M = asarray(M)
    if M.ndim != 2:
        raise ValueError("Wrong number of dimensions. It should be two.")
    M = DataFrame(M, columns=list(range(M.shape[1])))

    return M


def print_analysis(lik, name):
    lik_name = lik[0].upper() + lik[1:]
    print("*** %s using %s-GLMM ***" % (name, lik_name))


def phenotype_process(lik, y):
    if lik == 'poisson':
        y = clip(y, 0., 25000.)

    if isinstance(y, (tuple, list)):
        y = tuple([asarray(p, float) for p in y])
    else:
        y = asarray(y, float)

    if not ddall(isfinite(y)):
        msg = "One or more values of the provided phenotype "
        msg += "is not finite."
        raise ValueError(msg)
    return y


def covariates_process(M, nsamples):
    from pandas import DataFrame

    if M is None:
        M = DataFrame({'offset': ones(nsamples, float)})
    else:
        M = assure_named(M)

    if not ddall(isfinite(M.values)):
        msg = "One or more values of the provided covariates "
        msg += "is not finite."
        raise ValueError(msg)

    return M


_cache = dict(K=dict(K=None, QS=None, hash=None))


def kinship_process(K, nsamples, verbose):

    if K is None:
        K = eye(nsamples)

    K = asarray(K)

    ah = array_hash(K)
    nvalid = _cache['K']['K'] is None or (_cache['K']['hash'] != ah)

    if nvalid:

        if not ddall(isfinite(K)):
            msg = "One or more values of the provided covariance matrix "
            msg += "is not finite."
            raise ValueError(msg)

        K = gower_norm(K)

        desc = "Eigen decomposition of the covariance matrix..."
        with Timer(desc=desc, disable=not verbose):
            QS = economic_qs(K)
            _cache['K']['hash'] = ah
            _cache['K']['QS'] = QS
            _cache['K']['K'] = K
    else:
        QS = _cache['K']['QS']
        K = _cache['K']['K']

    return K, QS
