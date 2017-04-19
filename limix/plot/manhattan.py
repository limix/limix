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

def plot_manhattan(posCum,pv,chromBounds=None,
                    thr=None,qv=None,lim=None,xticklabels=True,
                    alphaNS=0.1,alphaS=0.5,colorNS='DarkBlue',
                    colorS='Orange',ax=None,thr_plotting=None,
                    labelS=None,labelNS=None):
    r"""Produce a manhattan plot

    Args:
        posCum: cumulative position.
        pv (array_like): pvalues.
        chromBounds (array_like): chrom boundaries (optionally).
                                  If not supplied, everything will
                                  be plotted into a single chromosome.
        qv (array_like): qvalues.
                         If provided, threshold for significance is set
                         on qvalues but pvalues are plotted.
        thr (float): threshold for significance.
                     The default is 0.01 significance level
                     (bonferroni-adjusted) if qvs are not specified,
                     or 0.01 FDR if qvs specified.
        lim (float): top limit on y-axis.
                     The default value is -1.2*log(pv.min()).
        xticklabels (bool): if true, xtick labels are printed.
                            The default value is True.
        alphaNS (float): transparency value of non-significant variants.
                         Must be in [0, 1].
        alphaS (float): transparency of significant variants.
                        Must be in [0, 1].
        ax (:class:`matplotlib.axes.AxesSubplot`):
                the target handle for this figure.
                If None, the current axes is set.
        thr_plotting (float): if specified, only P-values that are smaller
                              than thr_plotting are plotted.
        labelS (str): optional plotting label for significant variants.
        labelNS (str): optional plotting label for non significnat loci.

    Returns:
        :class:`matplotlib.axes.AxesSubplot`: matplotlib subplot 

    Example
    -------

        .. plot::

            from numpy.random import RandomState
            from numpy import arange
            from limix.plot import plot_manhattan 
            from matplotlib import pyplot as plt
            random = RandomState(1)

            pv = random.rand(5000)
            pv[1200:1250] = random.rand(50)**4 
            posCum = arange(5000)
            chromBounds = arange(0, 5000, 1000)
            
            fig = plt.figure(1, figsize=(8,3))
            plt.subplot(111)
            plot_manhattan(posCum, pv, chromBounds=chromBounds) 
            plt.tight_layout()
            plt.show()
    """
    import matplotlib.pylab as plt
    if ax is None:
        ax = plt.gca()

    if thr==None:
        thr = 0.01/float(posCum.shape[0])

    if lim==None:
        lim=-1.2*sp.log10(sp.minimum(pv.min(),thr))

    if chromBounds is None:
        chromBounds = sp.array([[0,posCum.max()]])
    else:
        chromBounds = sp.concatenate([chromBounds,sp.array([posCum.max()])])

    n_chroms = chromBounds.shape[0]
    for chrom_i in range(0,n_chroms-1,2):
        plt.fill_between(posCum,0,lim,where=(posCum>chromBounds[chrom_i]) & (posCum<chromBounds[chrom_i+1]),facecolor='LightGray',linewidth=0,alpha=0.5)

    if thr_plotting is not None:
        if pv is not None:
            i_small = pv<thr_plotting
        elif qv is not None:
            i_small = qv<thr_plotting

        if qv is not None:
            qv = qv[i_small]
        if pv is not None:
            pv = pv[i_small]
        if posCum is not None:
            posCum=posCum[i_small]

    if qv==None:
        Isign = pv<thr
    else:
        Isign = qv<thr

    plt.plot(posCum[~Isign],-sp.log10(pv[~Isign]),'.',color=colorNS,ms=5,alpha=alphaNS,label=labelNS)
    plt.plot(posCum[Isign], -sp.log10(pv[Isign]), '.',color=colorS,ms=5,alpha=alphaS,label=labelS)

    if qv is not None:
        plt.plot([0,posCum.max()],[-sp.log10(thr),-sp.log10(thr)],'--',color='Gray')

    plt.ylim(0,lim)

    plt.ylabel('-log$_{10}$pv')
    plt.xlim(0,posCum.max())
    xticks = sp.array([chromBounds[i:i+2].mean() for i in range(chromBounds.shape[0]-1)])
    ax.set_xticks(xticks)
    plt.xticks(fontsize=6)

    if xticklabels:
        ax.set_xticklabels(sp.arange(1,n_chroms+1))
        plt.xlabel('Chromosome')
    else:
        ax.set_xticklabels([])

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    return ax

