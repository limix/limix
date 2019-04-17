import warnings

from .._data import conform_dataset, normalize_likelihood
from .._display import session_block


class VarDec(object):
    """
    Variance decompositon through GLMMs.

    Example
    -------

    .. doctest::

        >>> from limix.vardec import VarDec
        >>> from limix.stats import multivariate_normal as mvn
        >>> from numpy import ones, eye, concatenate, zeros, exp
        >>> from numpy.random import RandomState
        >>>
        >>> random = RandomState(0)
        >>> nsamples = 20
        >>>
        >>> M = random.randn(nsamples, 2)
        >>> M = (M - M.mean(0)) / M.std(0)
        >>> M = concatenate((ones((nsamples, 1)), M), axis=1)
        >>>
        >>> K0 = random.randn(nsamples, 10)
        >>> K0 = K0 @ K0.T
        >>> K0 /= K0.diagonal().mean()
        >>> K0 += eye(nsamples) * 1e-4
        >>>
        >>> K1 = random.randn(nsamples, 10)
        >>> K1 = K1 @ K1.T
        >>> K1 /= K1.diagonal().mean()
        >>> K1 += eye(nsamples) * 1e-4
        >>>
        >>> y = M @ random.randn(3) + mvn(random, zeros(nsamples), K0)
        >>> y += mvn(random, zeros(nsamples), K1)
        >>>
        >>> vardec = VarDec(y, "normal", M)
        >>> vardec.append(K0)
        >>> vardec.append(K1)
        >>> vardec.append_iid()
        >>>
        >>> vardec.fit(verbose=False)
        >>> print(vardec) # doctest: +FLOAT_CMP
        Variance decomposition
        ======================
        <BLANKLINE>
        𝐲 ~ 𝓝(𝙼𝜶, 0.385⋅𝙺 + 1.184⋅𝙺 + 0.000⋅𝙸)
        >>> y = exp((y - y.mean()) / y.std())
        >>> vardec = VarDec(y, "poisson", M)
        >>> vardec.append(K0)
        >>> vardec.append(K1)
        >>> vardec.append_iid()
        >>>
        >>> vardec.fit(verbose=False)
        >>> print(vardec) # doctest: +FLOAT_CMP
        Variance decomposition
        ======================
        <BLANKLINE>
        𝐳 ~ 𝓝(𝙼𝜶, 0.000⋅𝙺 + 0.350⋅𝙺 + 0.000⋅𝙸) for yᵢ ~ Poisson(λᵢ=g(zᵢ)) and g(x)=eˣ
    """

    def __init__(self, y, lik="normal", M=None):
        """
        Constructor.

        Parameters
        ----------
        y : array_like
            Phenotype.
        lik : tuple, "normal", "bernoulli", "probit", "binomial", "poisson"
            Sample likelihood describing the residual distribution.
            Either a tuple or a string specifying the likelihood is required. The
            Normal, Bernoulli, Probit, and Poisson likelihoods can be selected by
            providing a string. Binomial likelihood on the other hand requires a tuple
            because of the number of trials: ``("binomial", array_like)``. Defaults to
            ``"normal"``.
        M : n×c array_like
            Covariates matrix.
        """
        from numpy import asarray
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
        """
        Covariace effect sizes.

        Returns
        -------
        effsizes : ndarray
            Effect sizes.
        """
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
        with session_block("Variance decomposition", disable=not verbose):
            if self._lik[0] == "normal":
                if self._simple_model():
                    self._fit_lmm_simple_model(verbose)
                else:
                    self._fit_lmm(verbose)
            else:
                if self._simple_model():
                    self._fit_glmm_simple_model(verbose)
                else:
                    self._fit_glmm(verbose)

            if verbose:
                print(self)

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
        from numpy_sugar import is_all_finite
        from numpy import asarray
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

    def plot(self):
        import limix
        import seaborn as sns
        from matplotlib.ticker import FormatStrFormatter

        variances = [c.scale for c in self._covariance]
        variances = [(v / sum(variances)) * 100 for v in variances]
        names = [c.name for c in self._covariance]

        ax = sns.barplot(x=names, y=variances)
        ax.yaxis.set_major_formatter(FormatStrFormatter("%.0f%%"))
        ax.set_xlabel("random effects")
        ax.set_ylabel("explained variance")
        ax.set_title("Variance decomposition")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            limix.plot.get_pyplot().tight_layout()
        limix.plot.show()

    def _fit_lmm(self, verbose):
        from glimix_core.cov import SumCov
        from glimix_core.gp import GP
        from numpy import asarray

        y = asarray(self._y, float).ravel()
        gp = GP(y, self._mean, SumCov(self._covariance))
        gp.fit(verbose=verbose)
        self._glmm = gp

    def _fit_glmm(self, verbose):
        from glimix_core.cov import SumCov
        from glimix_core.ggp import ExpFamGP
        from numpy import asarray

        y = asarray(self._y, float).ravel()
        gp = ExpFamGP(y, self._lik, self._mean, SumCov(self._covariance))
        gp.fit(verbose=verbose)
        self._glmm = gp

    def _fit_lmm_simple_model(self, verbose):
        from numpy_sugar.linalg import economic_qs
        from glimix_core.lmm import LMM
        from numpy import asarray

        K = self._get_matrix_simple_model()

        y = asarray(self._y, float).ravel()
        QS = None
        if K is not None:
            QS = economic_qs(K)
        lmm = LMM(y, self._M, QS)
        lmm.fit(verbose=verbose)
        self._set_simple_model_variances(lmm.v0, lmm.v1)
        self._glmm = lmm

    def _fit_glmm_simple_model(self, verbose):
        from numpy_sugar.linalg import economic_qs
        from glimix_core.glmm import GLMMExpFam
        from numpy import asarray

        K = self._get_matrix_simple_model()

        y = asarray(self._y, float).ravel()
        QS = None
        if K is not None:
            QS = economic_qs(K)

        glmm = GLMMExpFam(y, self._lik, self._M, QS)
        glmm.fit(verbose=verbose)

        self._set_simple_model_variances(glmm.v0, glmm.v1)
        self._glmm = glmm

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

    def __repr__(self):
        from glimix_core.cov import GivenCov
        from limix.qtl._result._draw import draw_model, draw_title

        covariance = ""
        for c in self._covariance:
            s = c.scale
            if isinstance(c, GivenCov):
                covariance += f"{s:.3f}⋅𝙺 + "
            else:
                covariance += f"{s:.3f}⋅𝙸 + "
        if len(covariance) > 2:
            covariance = covariance[:-3]

        msg = draw_title("Variance decomposition")
        msg += draw_model(self._lik[0], "𝙼𝜶", covariance)
        msg = msg.rstrip()

        return msg
