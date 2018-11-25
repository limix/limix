import sys

from . import _cov as user_cov, _mean as user_mean
from .. import _display
from .._data import check_likelihood_name
from ..qc._lik import normalise_extreme_values


if sys.version_info < (3, 0):
    PY2 = True
else:
    PY2 = False


class GLMMComposer(object):
    """ Construct GLMMs with any number of fixed and random effects.
    """

    def __init__(self, nsamples):
        """ Build a stub GLMM with a given number of samples.

        Parameters
        ----------
        nsamples : int
            Number of samples.
        """
        self._nsamples = nsamples
        self._likname = "normal"
        self._y = None
        self._fixed_effects = FixedEffects(nsamples)
        self._covariance_matrices = CovarianceMatrices(nsamples)
        self._glmm = None

    def decomp(self):
        r""" Get the fixed and random effects.

        Returns
        -------
        fixed_effects : Fixed effects.
        random_effects : Random effects.
        """
        from numpy import var as npvar

        decomp = dict(fixed_effects={}, random_effects={})

        for fe in self.fixed_effects:
            if hasattr(fe, "offset"):
                continue
            decomp["fixed_effects"][fe.name] = npvar(fe.value())

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
        """ Get likelihood name.

        Returns
        -------
        str
            Likelihood name.
        """
        return self._likname

    @likname.setter
    def likname(self, likname):
        """ Set likelihood name.

        Parameters
        ----------
        likname : str
            Likelihood name.
        """
        check_likelihood_name(likname)
        self._likname = likname.lower()
        self._glmm = None

    @property
    def y(self):
        """ Get the outcome array.

        Returns
        -------
        DataArray
            Outcome array.
        """
        return self._y

    @y.setter
    def y(self, y):
        """ Set the outcome array.

        Parameters
        ----------
        y : array_like
            Outcome array.
        """
        from numpy import all as npall, isfinite

        if not npall(isfinite(y)):
            raise ValueError("Phenotype values must be finite.")
        self._glmm = None
        self._y = normalise_extreme_values(y, "normal")

    @property
    def fixed_effects(self):
        """ Get the fixed effects.

        Returns
        -------
        dict
            Fixed effects.
        """
        return self._fixed_effects

    @property
    def covariance_matrices(self):
        """ Get the covariance matrices.

        Returns
        -------
        dict
            Covariance matrices.
        """
        return self._covariance_matrices

    def fit(self, verbose=True):
        """ Fit the model.

        Parameters
        ----------
        verbose : bool, optional
            Set ``False`` to silence it. Defaults to ``True``.
        """
        if self._likname == "normal":
            session_name = "composed lmm"
        else:
            session_name = "composed {}-glmm".format(self._likname)
        with _display.session_block(session_name, disable=not verbose):
            self._build_glmm()
            self._glmm.fit(verbose=verbose)

            if verbose:
                sys.stdout.flush()
                txt = _display.bold(str(self))
                _display.display(_display.format_richtext(txt))

    def lml(self):
        """ Get the log of the marginal likelihood.

        Returns
        -------
        float
            Log of the marginal likelihood.
        """
        self._build_glmm()
        return self._glmm.lml()

    def _build_glmm(self):
        from glimix_core.gp import GP
        from numpy import asarray

        if self._y is None:
            raise ValueError("Phenotype has not been set.")

        if self._likname == "normal" and self._glmm is None:
            gp = GP(
                asarray(self._y, float).ravel(),
                self._fixed_effects.impl,
                self._covariance_matrices.impl,
            )
            self._glmm = gp
            return

        if self._likname != "normal":
            raise NotImplementedError()

    def __repr__(self):
        from textwrap import TextWrapper

        width = _display.width()

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
        from numpy import arange

        self._sample_idx = arange(nsamples)
        self._fixed_effects = {"impl": [], "user": []}
        self._mean = None

    def __len__(self):
        return len(self._fixed_effects["impl"])

    def __getitem__(self, i):
        return self._fixed_effects["user"][i]

    def _setup_mean(self):
        from glimix_core.mean import SumMean

        if self._mean is None:
            mean = SumMean(self._fixed_effects["impl"])
            self._mean = {"impl": mean, "user": user_mean.SumMean(mean)}

    @property
    def impl(self):
        self._setup_mean()
        return self._mean["impl"]

    def append_offset(self):
        from glimix_core.mean import OffsetMean

        mean = OffsetMean()
        mean.set_data(self._sample_idx)
        self._fixed_effects["impl"].append(mean)
        self._fixed_effects["user"].append(user_mean.OffsetMean(mean))
        self._fixed_effects["user"][-1].name = "offset"
        self._mean = None

    def append(self, m, name=None):
        from numpy import all as npall, asarray, atleast_2d, isfinite
        from glimix_core.mean import LinearMean

        m = asarray(m, float)
        if m.ndim > 2:
            raise ValueError("Fixed-effect has to have between one and two dimensions.")

        if not npall(isfinite(m)):
            raise ValueError("Fixed-effect values must be finite.")

        m = atleast_2d(m.T).T
        mean = LinearMean(m.shape[1])
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
        from numpy import asarray

        vals = []
        for fi in self._fixed_effects["user"]:
            if isinstance(fi, user_mean.OffsetMean):
                vals.append(fi.offset)
            else:
                vals += list(fi.effsizes)
        return str(asarray(vals, float))


class CovarianceMatrices(object):
    def __init__(self, nsamples):
        from numpy import arange

        self._sample_idx = arange(nsamples)
        self._covariance_matrices = {"impl": [], "user": []}
        self._cov = None

    def __len__(self):
        return len(self._covariance_matrices["impl"])

    def __getitem__(self, i):
        return self._covariance_matrices["user"][i]

    def _setup_cov(self):
        from glimix_core.cov import SumCov

        if self._cov is None:
            cov = SumCov(self._covariance_matrices["impl"])
            self._cov = {"impl": cov, "user": user_cov.SumCov(cov)}

    @property
    def impl(self):
        self._setup_cov()
        return self._cov["impl"]

    def append_iid_noise(self):
        from glimix_core.cov import EyeCov

        cov = EyeCov()
        cov.set_data((self._sample_idx, self._sample_idx))
        self._covariance_matrices["impl"].append(cov)
        self._covariance_matrices["user"].append(user_cov.EyeCov(cov))
        self._covariance_matrices["user"][-1].name = "residual"
        self._cov = None

    def append(self, K, name=None):
        from numpy import all as npall, isfinite, issubdtype, number
        from glimix_core.cov import GivenCov

        if not issubdtype(K.dtype, number):
            raise ValueError("covariance-matrix is not numeric.")

        if K.ndim != 2:
            raise ValueError("Covariance-matrix has to have two dimensions.")

        if not npall(isfinite(K)):
            raise ValueError("Covariance-matrix values must be finite.")

        cov = GivenCov(K)
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
        from numpy import asarray

        vals = []
        for cm in self._covariance_matrices["user"]:
            vals.append(cm.scale)
        return str(asarray(vals, float))
