from __future__ import division

import bokeh
import bokeh.plotting
from numpy import (append, arange, argsort, flipud, inf, linspace, log10,
                   logspace, partition, searchsorted)
from scipy.special import betaincinv

def qqplot2(df,
           figure=None,
           colors=None,
           show=True,
           tools=None,
           nmax_points=1000,
           atleast_points=0.01,
           significance_level=0.01,
           paper_settings=False,
           **kwargs):
    r"""Plot number of significant hits across p-value thresholds.

    Args:
        df (:class:`pandas.DataFrame`): Columns `label` and `p-value`
            define labeled curves.

    Example
    -------

        .. bokeh-plot::

            from limix.plot import qqplot

            import pandas as pd
            import numpy as np
            random = np.random.RandomState(0)

            snp_ids = np.arange(1000)

            data1 = np.stack((['method1']*1000, random.rand(1000)),
                             axis=1)
            df1 = pd.DataFrame(data1, columns=['label', 'p-value'],
                               index=snp_ids)

            data2 = np.stack((['method2']*1000, random.rand(1000)),
                             axis=1)
            df2 = pd.DataFrame(data2, columns=['label', 'p-value'],
                               index=snp_ids)

            df = pd.concat([df1, df2])

            qqplot(df)
    """

    assert nmax_points > 1


    if tools is None:
        tools = ['save']

    if figure is None:
        figure = bokeh.plotting.figure(
            title='qqplot',
            tools=tools,
            x_axis_label="theoretical -log10(p-value)",
            y_axis_label="observed -log10(p-value)",
            **kwargs)

    labels = df['label'].unique()
    colors = _get_colors(colors, labels)
    threshold = _threshold(labels, df, nmax_points, atleast_points)

    npvals = inf

    for label in labels:
        pv = df[df['label'] == label]['p-value'].astype(float).values
        pv.sort()

        npvals = min(npvals, len(pv))

        lpv = -log10(flipud(pv))
        expected_lpv = _expected(len(lpv))

        i = searchsorted(pv, threshold)

        figure.circle(
            expected_lpv[-i:],
            lpv[-i:],
            color=colors[label],
            line_width=0,
            line_color=None,
            legend=label)

    _plot_confidence_band(npvals, nmax_points, atleast_points, figure,
                          significance_level)

    if paper_settings:
        _set_figure_for_paper(figure)

    if show:
        bokeh.plotting.show(figure)

    return figure


def _expected(n):
    lnpv = linspace(1 / (n + 1), n / (n + 1), n, endpoint=True)
    return flipud(-log10(lnpv))


def _rank_confidence_band(nranks, nmax_points, atleast_points,
                          significance_level):
    alpha = significance_level

    if nmax_points > nranks - 1:
        k0 = arange(1, nranks + 1)
    else:
        npoints = max(nmax_points, int(atleast_points * (nranks + 1)))
        k0 = linspace(1, nranks + 1, nmax_points, dtype=int)

    k1 = flipud(k0).copy()

    top = betaincinv(k0, k1, 1 - alpha)
    mean = k0 / (nranks + 1)
    bottom = betaincinv(k0, k1, alpha)

    return (bottom, mean, top)


def _labels(df):
    level = df.index.names.index('label')
    assert level == 0
    return list(df.index.get_level_values(level).unique())


def _threshold(labels, df, nmax_points, atleast_points):
    thr = 1.0
    for label in labels:
        pv = df[df['label'] == label]['p-value'].astype(float).values
        pv.sort()
        if len(pv) > nmax_points:
            npoints = max(nmax_points, int(atleast_points * len(pv)))
            npoints = min(len(pv) - 1, npoints)
            thr = min(thr, (pv[npoints - 1] + pv[npoints]) / 2)
    return thr


def _plot_confidence_band(npvals, nmax_points, atleast_points, figure,
                          significance_level):

    (bo, me, to) = _rank_confidence_band(npvals, nmax_points, atleast_points,
                                         significance_level)

    bo = flipud(-log10(bo))
    me = flipud(-log10(me))
    to = flipud(-log10(to))

    figure.line([me[-1], me[0]], [me[-1], me[0]], color='black')
    figure.legend.location = 'top_left'

    band_x = append(me, me[::-1])
    band_y = append(bo, to[::-1])
    figure.patch(
        band_x,
        band_y,
        line_color='black',
        fill_color='black',
        fill_alpha=0.15,
        line_alpha=0.5)

def _set_figure_for_paper(figure):
    figure.xaxis.axis_label_text_font_size = "24pt"
    figure.yaxis.axis_label_text_font_size = "24pt"
    figure.legend.label_text_font_size = "22pt"
    figure.xaxis.major_label_text_font_size = "18pt"
    figure.yaxis.major_label_text_font_size = "18pt"

def _get_colors(colors, labels):
    from bokeh.palettes import brewer

    if colors is None:
        colors = dict()

        colors_iter = iter(brewer['Spectral'][11])
        for label in labels:
            colors[label] = next(colors_iter)
    return colors

def _qqplot_bar(M=1000000, alphaLevel = 0.05, distr = 'log10'):
	"""calculate theoretical expectations for qqplot"""
	mRange=10**(sp.arange(sp.log10(0.5),sp.log10(M-0.5)+0.1,0.1));#should be exp or 10**?
	numPts=len(mRange);
	betaalphaLevel=sp.zeros(numPts);#down in the plot
	betaOneMinusalphaLevel=sp.zeros(numPts);#up in the plot
	betaInvHalf=sp.zeros(numPts);
	for n in range(numPts):
	   m=mRange[n]; #numPLessThanThresh=m;
	   betaInvHalf[n]=st.beta.ppf(0.5,m,M-m);
	   betaalphaLevel[n]=st.beta.ppf(alphaLevel,m,M-m);
	   betaOneMinusalphaLevel[n]=st.beta.ppf(1-alphaLevel,m,M-m);
	betaDown=betaInvHalf-betaalphaLevel;
	betaUp=betaOneMinusalphaLevel-betaInvHalf;

	theoreticalPvals=mRange/M;
	return betaUp, betaDown, theoreticalPvals


def qqplot(pv, distr = 'log10', alphaLevel = 0.05):
	"""
	This script makes a Quantile-Quantile plot of the observed
	negative log P-value distribution against the theoretical one under the null.

	Input:
		pv				pvalues (numpy array)
		distr           scale of the distribution (log10 or chi2)
		alphaLevel      significance bounds
	"""
	import matplotlib.pylab as plt
	shape_ok = (len(pv.shape)==1) or ((len(pv.shape)==2) and pv.shape[1]==1)
	assert shape_ok, 'qqplot requires a 1D array of p-values'

	tests = pv.shape[0]
	pnull = (0.5 + sp.arange(tests))/tests
	# pnull = np.sort(np.random.uniform(size = tests))
	Ipv = sp.argsort(pv)

	if distr == 'chi2':
	    qnull = sp.stats.chi2.isf(pnull, 1)
	    qemp = (sp.stats.chi2.isf(pv[Ipv],1))
	    xl = 'LOD scores'
	    yl = '$\chi^2$ quantiles'

	if distr == 'log10':
	    qnull = -sp.log10(pnull)
	    qemp = -sp.log10(pv[Ipv])

	    xl = '-log10(P) observed'
	    yl = '-log10(P) expected'

	plt.plot(qnull, qemp, '.')
	#plt.plot([0,qemp.m0x()], [0,qemp.max()],'r')
	plt.plot([0,qnull.max()], [0,qnull.max()],'r')
	plt.ylabel(xl)
	plt.xlabel(yl)
	if alphaLevel is not None:
	    if distr == 'log10':
	        betaUp, betaDown, theoreticalPvals = _qqplot_bar(M=tests,alphaLevel=alphaLevel,distr=distr)
	        lower = -sp.log10(theoreticalPvals-betaDown)
	        upper = -sp.log10(theoreticalPvals+betaUp)
	        plt.fill_between(-sp.log10(theoreticalPvals),lower,upper,color='grey',alpha=0.5)
	        #plt.plot(-sp.log10(theoreticalPvals),lower,'g-.')
	        #plt.plot(-sp.log10(theoreticalPvals),upper,'g-.')
