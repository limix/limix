import scipy as sp
import scipy.stats as st


def _qqplot_bar(M=1000000, alphaLevel=0.05, distr='log10'):
    """calculate theoretical expectations for qqplot"""
    mRange = 10**(sp.arange(sp.log10(0.5), sp.log10(M - 0.5) +
                            0.1, 0.1))  # should be exp or 10**?
    numPts = len(mRange)
    betaalphaLevel = sp.zeros(numPts)  # down in the plot
    betaOneMinusalphaLevel = sp.zeros(numPts)  # up in the plot
    betaInvHalf = sp.zeros(numPts)
    for n in range(numPts):
        m = mRange[n]  # numPLessThanThresh=m
        betaInvHalf[n] = st.beta.ppf(0.5, m, M - m)
        betaalphaLevel[n] = st.beta.ppf(alphaLevel, m, M - m)
        betaOneMinusalphaLevel[n] = st.beta.ppf(1 - alphaLevel, m, M - m)
    betaDown = betaInvHalf - betaalphaLevel
    betaUp = betaOneMinusalphaLevel - betaInvHalf

    theoreticalPvals = mRange / M
    return betaUp, betaDown, theoreticalPvals


def qqplot(pv, label='unknown', distr='log10', alphaLevel=0.05, ax=None,
           color=None):
    r"""Produces a Quantile-Quantile plot of the observed P value
        distribution against the theoretical one under the null.

    Parameters
    ----------

    pv : array_like
        P-values.
    distr : {'log10', 'chi2'}
        Scale of the distribution. If 'log10' is specified, the distribution
        of the -log10 P values is considered.
        If the distribution of the corresponding chi2-distributed test
        statistics is considered. Defaults to 'log10'.
    alphaLevel : float
        Significance bound.
    ax : :class:`matplotlib.axes.AxesSubplot`
        The target handle for this figure. If None, the current axes is set.

    Returns
    -------
    :class:`matplotlib.axes.AxesSubplot`
        Axes.

    Examples
    --------

    .. plot::

        from limix.plot import qqplot
        from numpy.random import RandomState
        from matplotlib import pyplot as plt
        random = RandomState(1)

        pv = random.rand(10000)

        fig = plt.figure(1, figsize=(5,5))
        plt.subplot(111)
        qqplot(pv)
        plt.tight_layout()
        plt.show()
    """
    import matplotlib.pylab as plt
    if ax is None:
        ax = plt.gca()

    shape_ok = (len(pv.shape) == 1) or (
        (len(pv.shape) == 2) and pv.shape[1] == 1)
    assert shape_ok, 'qqplot requires a 1D array of p-values'

    tests = pv.shape[0]
    pnull = (0.5 + sp.arange(tests)) / tests
    # pnull = np.sort(np.random.uniform(size = tests))
    Ipv = sp.argsort(pv)

    if distr == 'chi2':
        qnull = sp.stats.chi2.isf(pnull, 1)
        qemp = (sp.stats.chi2.isf(pv[Ipv], 1))
        xl = 'LOD scores'
        yl = '$\chi^2$ quantiles'

    if distr == 'log10':
        qnull = -sp.log10(pnull)
        qemp = -sp.log10(pv[Ipv])

        xl = '-log10(P) observed'
        yl = '-log10(P) expected'

    plt.plot(qnull, qemp, '.', color=color, label=label)
    # plt.plot([0,qemp.m0x()], [0,qemp.max()],'r')
    plt.plot([0, qnull.max()], [0, qnull.max()], 'r')
    plt.ylabel(xl)
    plt.xlabel(yl)
    if alphaLevel is not None:
        if distr == 'log10':
            betaUp, betaDown, theoreticalPvals = _qqplot_bar(
                M=tests, alphaLevel=alphaLevel, distr=distr)
            lower = -sp.log10(theoreticalPvals - betaDown)
            upper = -sp.log10(theoreticalPvals + betaUp)
            plt.fill_between(-sp.log10(theoreticalPvals),
                             lower, upper, color='grey', alpha=0.5)
            # plt.plot(-sp.log10(theoreticalPvals),lower,'g-.')
            # plt.plot(-sp.log10(theoreticalPvals),upper,'g-.')

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    return ax
