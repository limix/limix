from numpy import all as npall
from numpy import arange, asarray, atleast_2d, isfinite
from textwrap import TextWrapper

import glimix_core
from glimix_core.gp import GP
from blessings import Terminal

from . import _cov as user_cov
from . import _mean as user_mean
from ..likelihood import assert_likelihood_name
from ..util import session_text


class GLMMComposer(object):
    def __init__(self, nsamples):
        self._nsamples = nsamples
        self._likname = "normal"
        self._y = None
        self._fixed_effects = FixedEffects(nsamples)
        self._covariance_matrices = CovarianceMatrices(nsamples)
        self._glmm = None

    @property
    def likname(self):
        return self._likname

    @likname.setter
    def likname(self, likname):
        assert_likelihood_name(likname)
        self._likname = likname.lower()
        self._glmm = None

    @property
    def y(self):
        return self._y

    @y.setter
    def y(self, y):
        if not npall(isfinite(y)):
            raise ValueError("Phenotype values must be finite.")
        self._glmm = None
        self._y = y

    @property
    def fixed_effects(self):
        return self._fixed_effects

    @property
    def covariance_matrices(self):
        return self._covariance_matrices

    def fit(self, verbose=True, progress=True):
        if self._likname == "normal":
            session_name = "composed lmm"
        else:
            session_name = "composed {}-glmm".format(self._likname)
        with session_text(session_name):
            self._build_glmm()
            self._glmm.feed().maximize(verbose=progress)
            if verbose:
                print(self)

    def lml(self):
        self._build_glmm()
        return self._glmm.lml()

    def _build_glmm(self):
        if self._likname == "normal" and self._glmm is None:
            gp = GP(self._y, self._fixed_effects.impl, self._covariance_matrices.impl)
            self._glmm = gp
            return

        if self._likname != "normal":
            raise NotImplementedError()

    def __str__(self):
        term = Terminal()
        width = term.width
        if width is None:
            width = 79

        if self._likname == "normal":
            s = "GLMMComposer using LMM\n"
            s += "-----------------------\n"
        else:
            s = "unknown"

        s += "LML: {}\n".format(self.lml())

        w = TextWrapper(initial_indent="", subsequent_indent=" " * 21, width=width)
        s += w.fill("Fixed-effect sizes: " + str(self._fixed_effects)) + "\n"

        w = TextWrapper(initial_indent="", subsequent_indent=" " * 27, width=width)
        s += w.fill("Covariance-matrix scales: " + str(self._covariance_matrices))
        return term.bold(s)


class FixedEffects(object):
    def __init__(self, nsamples):
        self._sample_idx = arange(nsamples)
        self._fixed_effects = {"impl": [], "user": []}
        self._mean = None

    def __len__(self):
        return len(self._fixed_effects["impl"])

    def __getitem__(self, i):
        return self._fixed_effects["user"][i]

    def _setup_mean(self):
        if self._mean is None:
            mean = glimix_core.mean.SumMean(self._fixed_effects["impl"])
            self._mean = {"impl": mean, "user": user_mean.SumMean(mean)}

    @property
    def impl(self):
        self._setup_mean()
        return self._mean["impl"]

    def append_offset(self):
        mean = glimix_core.mean.OffsetMean()
        mean.set_data(self._sample_idx)
        self._fixed_effects["impl"].append(mean)
        self._fixed_effects["user"].append(user_mean.OffsetMean(mean))
        self._mean = None

    def append(self, m):
        m = asarray(m, float)
        if m.ndim > 2:
            raise ValueError("Fixed-effect has to have between one and two dimensions.")

        if not npall(isfinite(m)):
            raise ValueError("Fixed-effect values must be finite.")

        m = atleast_2d(m.T).T
        mean = glimix_core.mean.LinearMean(m.shape[1])
        mean.set_data(m)

        self._fixed_effects["impl"].append(mean)
        self._fixed_effects["user"].append(user_mean.LinearMean(mean))
        self._mean = None

    @property
    def mean(self):
        self._setup_mean()
        return self._mean["user"]

    def __str__(self):
        vals = []
        for fi in self._fixed_effects["user"]:
            if isinstance(fi, user_mean.OffsetMean):
                vals.append(fi.offset)
            else:
                vals += list(fi.effsizes)
        return str(asarray(vals, float))


class CovarianceMatrices(object):
    def __init__(self, nsamples):
        self._sample_idx = arange(nsamples)
        self._covariance_matrices = {"impl": [], "user": []}
        self._cov = None

    def __len__(self):
        return len(self._covariance_matrices["impl"])

    def __getitem__(self, i):
        return self._covariance_matrices["user"][i]

    def _setup_cov(self):
        if self._cov is None:
            cov = glimix_core.cov.SumCov(self._covariance_matrices["impl"])
            self._cov = {"impl": cov, "user": user_cov.SumCov(cov)}

    @property
    def impl(self):
        self._setup_cov()
        return self._cov["impl"]

    def append_iid_noise(self):
        cov = glimix_core.cov.EyeCov()
        cov.set_data((self._sample_idx, self._sample_idx))
        self._covariance_matrices["impl"].append(cov)
        self._covariance_matrices["user"].append(user_cov.EyeCov(cov))
        self._cov = None

    def append(self, K):
        if K.ndim != 2:
            raise ValueError("Covariance-matrix has to have two dimensions.")

        if not npall(isfinite(K)):
            raise ValueError("Covariance-matrix values must be finite.")

        cov = glimix_core.cov.GivenCov(K)
        cov.set_data((self._sample_idx, self._sample_idx))
        self._covariance_matrices["impl"].append(cov)
        self._covariance_matrices["user"].append(user_cov.GivenCov(cov))
        self._cov = None

    @property
    def cov(self):
        self._setup_cov()
        return self._cov["user"]

    def __str__(self):
        vals = []
        for cm in self._covariance_matrices["user"]:
            vals.append(cm.scale)
        return str(asarray(vals, float))

