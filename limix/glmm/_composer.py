from numpy import all as npall
from numpy import isfinite, atleast_2d, asarray, arange
from glimix_core.cov import EyeCov, GivenCov, SumCov
from glimix_core.mean import OffsetMean, LinearMean, SumMean
from glimix_core.gp import GP
from glimix_core.lmm import LMM
from ..likelihood import assert_likelihood_name


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

    def fit(self, verbose=True):
        self._build_glmm()
        self._glmm.feed().maximize(verbose=verbose)

    def lml(self):
        self._build_glmm()
        return glmm.lml()

    def _build_glmm(self):
        if self._likname == "normal" and self._glmm is not None:
            gp = GP(self._y, self._fixed_effects.mean, self._covariance_matrices.cov)
            self._glmm = gp
            return

        raise NotImplementedError()

    def __str__(self):
        s = analysis_welcome(self._likname, "GLMMComposer") + "\n"
        s += "Fixed-effect sizes: " + str(self._fixed_effects) + "\n"
        s += "Covariance-matrix scales: " + str(self._covariance_matrices)
        return s


class FixedEffects(object):
    def __init__(self, nsamples):
        self._sample_idx = arange(nsamples)
        self._fixed_effects = []
        self._mean = None

    def __len__(self):
        return len(self._fixed_effects)

    def __getitem__(self, i):
        return self._fixed_effects[i]

    def append_offset(self):
        mean = OffsetMean()
        mean.set_data(self._sample_idx)
        self._fixed_effects.append(mean)
        self._mean = None

    def append(self, m):
        m = asarray(m, float)
        if m.ndim > 2:
            raise ValueError("Fixed-effect has to have between one and two dimensions.")

        if not npall(isfinite(m)):
            raise ValueError("Fixed-effect values must be finite.")

        m = atleast_2d(m.T).T
        mean = LinearMean(m.shape[1])
        mean.set_data(m)

        self._fixed_effects.append(mean)
        self._mean = None

    @property
    def mean(self):
        if self._mean is None:
            self._mean = SumMean(self._fixed_effects)
        return self._mean

    def __str__(self):
        vals = []
        for fi in self._fixed_effects:
            if isinstance(fi, OffsetMean):
                vals.append(fi.offset)
            else:
                vals += list(fi.effsizes)
        return str(asarray(vals, float))


class CovarianceMatrices(object):
    def __init__(self, nsamples):
        self._sample_idx = arange(nsamples)
        self._covariance_matrices = []
        self._cov = None

    def __len__(self):
        return len(self._covariance_matrices)

    def __getitem__(self, i):
        return self._covariance_matrices[i]

    def append_iid_noise(self):
        cov = EyeCov()
        cov.set_data((self._sample_idx, self._sample_idx))
        self._covariance_matrices.append(cov)
        self._cov = None

    def append(self, K):
        if K.ndim != 2:
            raise ValueError("Covariance-matrix has to have two dimensions.")

        if not npall(isfinite(K)):
            raise ValueError("Covariance-matrix values must be finite.")

        cov = GivenCov(K)
        cov.set_data((self._sample_idx, self._sample_idx))
        self._covariance_matrices.append(cov)
        self._cov = None

    @property
    def cov(self):
        if self._cov is None:
            self._cov = SumCov(self._covariance_matrices)
        return self._cov

    def __str__(self):
        vals = []
        for cm in self._covariance_matrices:
            vals.append(cm.scale)
        return str(asarray(vals, float))


def analysis_welcome(lik, name):
    lik_name = lik[0].upper() + lik[1:]
    return f"*** {name} using {lik_name}-GLMM ***"
