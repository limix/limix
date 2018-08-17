import sys
from textwrap import TextWrapper

import glimix_core
from glimix_core.gp import GP
from numpy import (
    all as npall,
    arange,
    asarray,
    atleast_2d,
    isfinite,
    issubdtype,
    number,
    var as np_var,
)

from . import _cov as user_cov, _mean as user_mean
from .. import display
from .._likelihood import assert_likelihood_name, normalise_extreme_values

if sys.version_info < (3, 0):
    PY2 = True
else:
    PY2 = False


class GLMMComposer(object):
    def __init__(self, nsamples):
        self._nsamples = nsamples
        self._likname = "normal"
        self._y = None
        self._fixed_effects = FixedEffects(nsamples)
        self._covariance_matrices = CovarianceMatrices(nsamples)
        self._glmm = None

    def plot(self):
        from pandas import DataFrame
        from matplotlib.ticker import FormatStrFormatter
        import seaborn as sns

        fe_vars = []
        fe_names = []
        for fe in self.fixed_effects:
            if hasattr(fe, "offset"):
                continue
            fe_vars.append(np_var(fe.value()))
            fe_names.append(fe.name)

        scales = [cm.scale for cm in self.covariance_matrices]
        re_names = [re.name for re in self.covariance_matrices]

        fe = DataFrame([fe_vars], columns=fe_names)
        re = DataFrame([scales], columns=re_names)

        fe = fe.div(fe.sum(axis=1), axis=0).mean(axis=0)
        fe *= 100

        re = re.div(re.sum(axis=1), axis=0).mean(axis=0)
        re *= 100

        ax = sns.barplot(x=re.index, y=re.values)
        ax.yaxis.set_major_formatter(FormatStrFormatter("%.0f%%"))

        return ax

    def decomp(self):
        decomp = dict(fixed_effects={}, random_effects={})

        for fe in self.fixed_effects:
            if hasattr(fe, "offset"):
                continue
            decomp["fixed_effects"][fe.name] = np_var(fe.value())

        for re in self.covariance_matrices:
            decomp["random_effects"][re.name] = re.scale

        total = 0
        for _, v in iter(decomp.items()):
            for _, vi in iter(v.items()):
                total += vi

        for k0, v in iter(decomp.items()):
            for k1, vi in iter(v.items()):
                decomp[k0][k1] = vi / total

        return decomp

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
        self._y = normalise_extreme_values(y, "normal")

    @property
    def fixed_effects(self):
        return self._fixed_effects

    @property
    def covariance_matrices(self):
        return self._covariance_matrices

    def fit(self, verbose=True):
        if self._likname == "normal":
            session_name = "composed lmm"
        else:
            session_name = "composed {}-glmm".format(self._likname)
        with display.session_text(session_name, disable=not verbose):
            self._build_glmm()
            self._glmm.fit(verbose=verbose)

            if verbose:
                sys.stdout.flush()
                txt = display.bold(str(self))
                display.display(display.format_richtext(txt))

    def lml(self):
        self._build_glmm()
        return self._glmm.lml()

    def _build_glmm(self):
        if self._y is None:
            raise ValueError("Phenotype has not been set.")

        if self._likname == "normal" and self._glmm is None:
            gp = GP(self._y, self._fixed_effects.impl, self._covariance_matrices.impl)
            self._glmm = gp
            return

        if self._likname != "normal":
            raise NotImplementedError()

    def __repr__(self):
        width = display.width()

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
        return s

    def __str__(self):
        if PY2:
            return self.__unicode__().encode("utf-8")
        return self.__repr__()

    def __unicode__(self):
        return self.__repr__()


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
        self._fixed_effects["user"][-1].name = "offset"
        self._mean = None

    def append(self, m, name=None):
        m = asarray(m, float)
        if m.ndim > 2:
            raise ValueError("Fixed-effect has to have between one and two dimensions.")

        if not npall(isfinite(m)):
            raise ValueError("Fixed-effect values must be finite.")

        m = atleast_2d(m.T).T
        mean = glimix_core.mean.LinearMean(m.shape[1])
        mean.set_data(m)

        n = len(self._fixed_effects["impl"])
        if name is None:
            name = "unnamed-fe-{}".format(n)
        self._fixed_effects["impl"].append(mean)
        self._fixed_effects["user"].append(user_mean.LinearMean(mean))
        self._fixed_effects["user"][-1].name = name
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
        self._covariance_matrices["user"][-1].name = "residual"
        self._cov = None

    def append(self, K, name=None):
        if not issubdtype(K.dtype, number):
            raise ValueError("covariance-matrix is not numeric.")

        if K.ndim != 2:
            raise ValueError("Covariance-matrix has to have two dimensions.")

        if not npall(isfinite(K)):
            raise ValueError("Covariance-matrix values must be finite.")

        cov = glimix_core.cov.GivenCov(K)
        cov.set_data((self._sample_idx, self._sample_idx))

        n = len(self._covariance_matrices["impl"])
        if name is None:
            name = "unnamed-re-{}".format(n)

        self._covariance_matrices["impl"].append(cov)
        self._covariance_matrices["user"].append(user_cov.GivenCov(cov))
        self._covariance_matrices["user"][-1].name = name
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
