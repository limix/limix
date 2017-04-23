# Copyright(c) 2014, The LIMIX developers (Christoph Lippert, Paolo Francesco Casale, Oliver Stegle)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import scipy as sp


def plot_manhattan(pv, position=None, posCum=None, chromBounds=None,
                   thr=None, qv=None, lim=None, xticklabels=True,
                   alphaNS=0.1, alphaS=0.5, colorNS='DarkBlue',
                   colorS='Orange', ax=None, thr_plotting=None,
                   labelS=None, labelNS=None):
    r"""Produce a manhattan plot.

    Args:
        pv (array_like): pvalues.

        position (list or :class:`pandas.DataFrame`, optional):
            positions in chromosome/chromosomal basepair position format.
            It can be specified as

            - list ``[chrom, pos]`` where ``chrom`` and ``pos`` are
              *ndarray* with chromosome values and basepair positions;
            - pandas DataFrame of chromosome values (key='chrom') and
              basepair positions (key='pos').

            Alternatively, variant positions can be specified by
            setting directly the cumulative position (``posCum``).

        posCum (array_like, optional):
            cumulative position.
            It has effect only if ``position`` is not specified.
            By default, ``posCum`` is set to ``range(S)`` where
            `S` is the number of variants.

        chromBounds (array_like, optional):
            chrom boundaries on cumulative positions.
            It has effect only if ``position`` is not specified.
            If neither ``position`` nor ``chromBounds`` are
            specified, chomosome boundaries are not plotted.

        qv (array_like, optional):
            qvalues.
            If provided, threshold for significance is set
            on qvalues but pvalues are plotted.

        thr (float, optional):
            threshold for significance.
            The default is 0.01 significance level (bonferroni-adjusted)
            if qvs are not specified, or 0.01 FDR if qvs specified.

        lim (float, optional):
            top limit on y-axis.
            The default value is -1.2*log(pv.min()).

        xticklabels (bool, optional):
            if true, xtick labels are printed.
            The default value is True.

        alphaNS (float, optional):
            transparency value of non-significant variants.
            Must be in [0, 1].

        alphaS (float, optional):
            transparency of significant variants. Must be in [0, 1].

        ax (:class:`matplotlib.axes.AxesSubplot`):
            the target handle for this figure.
            If None, the current axes is set.

        thr_plotting (float):
            if specified, only P-values that are smaller
            than thr_plotting are plotted.

        labelS (str):
            optional plotting label for significant variants.

        labelNS (str):
            optional plotting label for non significnat loci.

    Returns:
        :class:`matplotlib.axes.AxesSubplot`: matplotlib subplot

    Examples
    --------

        .. plot::

            from numpy.random import RandomState
            from numpy import arange, ones, kron
            from limix.plot import plot_manhattan
            from matplotlib import pyplot as plt
            random = RandomState(1)

            pv = random.rand(5000)
            pv[1200:1250] = random.rand(50)**4

            chrom  = kron(arange(1,6), ones(1000))
            pos = kron(ones(5), arange(1,1001))

            fig = plt.figure(1, figsize=(8,3))
            plt.subplot(111)
            plot_manhattan(pv, [chrom, pos])
            plt.tight_layout()
            plt.show()
    """
    import matplotlib.pylab as plt
    from limix.util import estCumPos
    if ax is None:
        ax = plt.gca()

    if position is not None:
        posCum, chromBounds = estCumPos(position, return_chromstart=True)
    elif posCum is None:
        posCum = sp.arange(len(pv))

    if thr is None:
        thr = 0.01 / float(posCum.shape[0])

    if lim is None:
        lim = -1.2 * sp.log10(sp.minimum(pv.min(), thr))

    if chromBounds is None:
        chromBounds = sp.array([[0, posCum.max()]])
    else:
        chromBounds = sp.concatenate([chromBounds, sp.array([posCum.max()])])

    n_chroms = chromBounds.shape[0]
    for chrom_i in range(0, n_chroms - 1, 2):
        plt.fill_between(posCum, 0, lim, where=(posCum > chromBounds[chrom_i]) & (
            posCum < chromBounds[chrom_i + 1]), facecolor='LightGray', linewidth=0, alpha=0.5)

    if thr_plotting is not None:
        if pv is not None:
            i_small = pv < thr_plotting
        elif qv is not None:
            i_small = qv < thr_plotting

        if qv is not None:
            qv = qv[i_small]
        if pv is not None:
            pv = pv[i_small]
        if posCum is not None:
            posCum = posCum[i_small]

    if qv is None:
        Isign = pv < thr
    else:
        Isign = qv < thr

    plt.plot(posCum[~Isign], -sp.log10(pv[~Isign]), '.',
             color=colorNS, ms=5, alpha=alphaNS, label=labelNS)
    plt.plot(posCum[Isign], -sp.log10(pv[Isign]), '.',
             color=colorS, ms=5, alpha=alphaS, label=labelS)

    if qv is not None:
        plt.plot([0, posCum.max()], [-sp.log10(thr), -
                                     sp.log10(thr)], '--', color='Gray')

    plt.ylim(0, lim)

    plt.ylabel('-log$_{10}$pv')
    plt.xlim(0, posCum.max())
    xticks = sp.array([chromBounds[i:i + 2].mean()
                       for i in range(chromBounds.shape[0] - 1)])
    ax.set_xticks(xticks)
    plt.xticks(fontsize=6)

    if xticklabels:
        ax.set_xticklabels(sp.arange(1, n_chroms + 1))
        plt.xlabel('Chromosome')
    else:
        ax.set_xticklabels([])

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    return ax
