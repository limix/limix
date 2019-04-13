import sys

from numpy import asarray
from . import _cov as user_cov, _mean as user_mean
from .. import _display
from .._data import asarray as _asarray, conform_dataset, normalize_likelihood
from ..qc._lik import normalise_extreme_values


class VarDec(object):
    """
    Construct GLMMs with any number of fixed and random effects.
    """

    def __init__(self, y, lik, M=None):
        """
        Build a stub GLMM with a given number of samples.

        Parameters
        ----------
        nsamples : int
            Number of samples.
        """
        from glimix_core.mean import LinearMean

        y = asarray(y, float)
        data = conform_dataset(y, M)
        y = data["y"]
        M = data["M"]
        self._y = y
        self._M = M
        self._lik = normalize_likelihood(lik)
        self._mean = LinearMean(asarray(M, float))
        self._covariance = []
        self._glmm = None
        self._fit = False
        self._unnamed = 0

    @property
    def effsizes(self):
        if not self._fit:
            self.fit()
        return self._mean.effsizes

    @property
    def covariance(self):
        """
        Get the covariance matrices.

        Returns
        -------
        covariances : list
            Covariance matrices.
        """
        return self._covariance

    def fit(self, verbose=True):
        """
        Fit the model.

        Parameters
        ----------
        verbose : bool, optional
            Set ``False`` to silence it. Defaults to ``True``.
        """
        if self._lik[0] == "normal":
            session_name = "composed lmm"
        else:
            session_name = "composed {}-glmm".format(self._lik[0])
        with _display.session_block(session_name, disable=not verbose):
            if self._lik[0] == "normal":
                if self._simple_model():
                    self._fit_lmm_simple_model(verbose)
                else:
                    self._fit_lmm(verbose)

            # if verbose:
            #     sys.stdout.flush()
            #     txt = _display.bold(str(self))
            #     _display.display(_display.format_richtext(txt))

        self._fit = True

    def lml(self):
        """
        Get the log of the marginal likelihood.

        Returns
        -------
        float
            Log of the marginal likelihood.
        """
        if not self._fit:
            self._glmm.fit()
        return self._glmm.lml()

    def append_iid(self, name="residual"):
        from glimix_core.cov import EyeCov

        c = EyeCov(self._y.shape[0])
        c.name = name
        self._covariance.append(c)

    def append(self, K, name=None):
        from numpy import all as npall, isfinite, issubdtype, number
        from numpy_sugar import is_all_finite
        from glimix_core.cov import GivenCov

        data = conform_dataset(self._y, K=K)
        K = asarray(data["K"], float)

        if not is_all_finite(K):
            raise ValueError("Covariance-matrix values must be finite.")

        K = K / K.diagonal().mean()
        cov = GivenCov(K)
        if name is None:
            name = "unnamed-{}".format(self._unnamed)
            self._unnamed += 1
        cov.name = name

        self._covariance.append(cov)

    def _fit_lmm(self, verbose):
        from glimix_core.cov import SumCov
        from glimix_core.gp import GP

        y = asarray(self._y, float).ravel()
        gp = GP(y, self._mean, SumCov(self._covariance))
        gp.fit(verbose=verbose)
        self._glmm = gp

    def _fit_lmm_simple_model(self, verbose):
        from numpy_sugar.linalg import economic_qs
        from glimix_core.lmm import LMM

        K = self._get_matrix_simple_model()

        y = asarray(self._y, float).ravel()
        QS = None
        if K is not None:
            QS = economic_qs(K)
        lmm = LMM(y, self._M, QS)
        lmm.fit(verbose=verbose)
        self._set_simple_model_variances(lmm.v0, lmm.v1)
        self._glmm = lmm

    def _set_simple_model_variances(self, v0, v1):
        from glimix_core.cov import GivenCov, EyeCov

        for c in self._covariance:
            if isinstance(c, GivenCov):
                c.scale = v0
            elif isinstance(c, EyeCov):
                c.scale = v1

    def _get_matrix_simple_model(self):
        from glimix_core.cov import GivenCov

        K = None
        for i in range(len(self._covariance)):
            if isinstance(self._covariance[i], GivenCov):
                self._covariance[i].scale = 1.0
                K = self._covariance[i].value()
                break
        return K

    def _fit_glmm(self, verbose):
        from glimix_core.gp import GP

    def _simple_model(self):
        from glimix_core.cov import GivenCov, EyeCov

        if len(self._covariance) > 2:
            return False

        c = self._covariance
        if len(c) == 1 and isinstance(c[0], EyeCov):
            return True

        if isinstance(c[0], GivenCov) and isinstance(c[1], EyeCov):
            return True

        if isinstance(c[1], GivenCov) and isinstance(c[0], EyeCov):
            return True

        return False

    # def _build_glmm(self):
    #     from numpy import asarray

    #     if self._y is None:
    #         raise ValueError("Phenotype has not been set.")

    #     if self._likname == "normal" and self._glmm is None:
    #         gp = GP(
    #             asarray(self._y, float).ravel(),
    #             self._fixed_effects.impl,
    #             self._covariance.impl,
    #         )
    #         self._glmm = gp
    #         return

    #     if self._likname != "normal":
    #         raise NotImplementedError()

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
        s += w.fill("Covariance-matrix scales: " + str(self._covariance))
        return s
