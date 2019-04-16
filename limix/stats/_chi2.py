class Chi2Mixture(object):
    """
    Mixture of 𝜒² distributions.

    Class for evaluation of P values for a test statistic that follows a
    two-component mixture of chi2

    .. math::

        p(x) = (1-𝑝)𝜒²(0) + 𝑝𝑎𝜒²(𝑑).

    Here 𝑝 is the probability being in the first component and 𝑎 and 𝑑 are the scale
    parameter and the number of degrees of freedom of the second
    component.

    .. warning::
        This class is not production-ready. Keep in mind that this interface is likely
        to change.

    Parameters
    ----------
    scale_min : float
        Minimum value used for fitting the scale parameter.
    scale_max : float
        Maximum value used for fitting the scale parameter.
    dofmin : float
        Minimum value used for fitting the dof parameter.
    dofmax : float
        Maximum value used for fitting the dof parameter.
    qmax : float
        Only the top qmax quantile is used for the fit.
    n_interval : int
        Number of intervals when performing gridsearch.
    tol : float
        Tolerance of being zero.

    Example
    -------

    .. doctest::

        >>> from numpy.random import RandomState
        >>> import scipy as sp
        >>> from limix.stats import Chi2Mixture
        >>>
        >>> scale = 0.3
        >>> dof = 2
        >>> mixture = 0.2
        >>> n = 100
        >>>
        >>> random = RandomState(1)
        >>> x =  random.chisquare(dof, n)
        >>> n0 = int( (1-mixture) * n)
        >>> idxs = random.choice(n, n0, replace=False)
        >>> x[idxs] = 0
        >>>
        >>> chi2mix = Chi2Mixture(scale_min=0.1, scale_max=5.0,
        ...                       dof_min=0.1, dof_max=5.0,
        ...                       qmax=0.1, tol=4e-3)
        >>> chi2mix.estimate_chi2mixture(x)
        >>> pv = chi2mix.sf(x)
        >>> print(pv[:4]) # doctest: +FLOAT_CMP
        [0.2 0.2 0.2 0.2]
        >>>
        >>> print('%.2f' % chi2mix.scale)
        1.98
        >>> print('%.2f' % chi2mix.dof)
        0.89
        >>> print('%.2f' % chi2mix.mixture)
        0.20
    """

    def __init__(
        self,
        scale_min=0.1,
        scale_max=5.0,
        dof_min=0.1,
        dof_max=5.0,
        n_intervals=100,
        qmax=0.1,
        tol=0,
        lrt=None,
    ):

        self.scale_min = scale_min
        self.scale_max = scale_max
        self.dof_min = dof_min
        self.dof_max = dof_max
        self.qmax = qmax
        self.n_intervals = n_intervals
        self.tol = tol

    def estimate_chi2mixture(self, lrt):
        r"""Estimates the parameters of a chi2 mixture.

        Estimates the parameters of a chi2 mixture by fitting the empirical
        distribution of null test statistic.

        Parameters
        ----------
        lrt : array_like
            Null test statistcs.
        """
        import scipy as sp
        import scipy.stats as st

        # step 1: estimate the probability of being in component one
        self.mixture = 1 - (lrt <= self.tol).mean()
        n_false = sp.sum(lrt > self.tol)

        # step 2: only use the largest qmax fraction of test statistics to
        #         estimate the remaining parameters
        n_fitting = int(sp.ceil(self.qmax * n_false))
        lrt_sorted = -sp.sort(-lrt)[:n_fitting]
        q = sp.linspace(0, 1, n_false)[1 : n_fitting + 1]
        log_q = sp.log10(q)

        # step 3: fitting scale and dof by minimizing the squared error
        #          of the log10 p-values with their theorietical values
        #          [uniform distribution]
        MSE_opt = sp.inf
        MSE = sp.zeros((self.n_intervals, self.n_intervals))

        for i, scale in enumerate(
            sp.linspace(self.scale_min, self.scale_max, self.n_intervals)
        ):
            for j, dof in enumerate(
                sp.linspace(self.dof_min, self.dof_max, self.n_intervals)
            ):
                p = st.chi2.sf(lrt_sorted / scale, dof)
                log_p = sp.log10(p)
                MSE[i, j] = sp.mean((log_q - log_p) ** 2)
                if MSE[i, j] < MSE_opt:
                    MSE_opt = MSE[i, j]
                    self.scale = scale
                    self.dof = dof

    def sf(self, lrt):
        r"""Computes the p-values from test statistics lrt.

        Parameters
        ----------
        lrt : array_like
            Test statistics.

        Returns
        -------
        array_like
            P-values.
        """
        from numpy import asarray
        import scipy as sp
        import scipy.stats as st

        lrt = asarray(lrt, float)

        _lrt = sp.copy(lrt)
        _lrt[lrt < self.tol] = 0
        pv = self.mixture * st.chi2.sf(_lrt / self.scale, self.dof)
        return pv
