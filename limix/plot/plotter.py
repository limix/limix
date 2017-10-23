import matplotlib as mpl
from matplotlib import pyplot as plt

_figure_id = 98748


def get_plot():
    plt.close(_figure_id)

    _set_rc_params()
    figure = plt.figure(_figure_id)
    axes = figure.add_subplot(111)

    return LimixPlot(figure, axes)


class LimixPlot(object):
    def __init__(self, figure, axes):
        self._figure = figure
        self._axes = axes
        axes.set_axisbelow(True)
        axes.spines['right'].set_visible(True)
        axes.spines['top'].set_visible(True)

    @property
    def figure(self):
        return self._figure

    @property
    def axes(self):
        return self._axes

    def ion(self):
        plt.ion()

    def ioff(self):
        plt.ioff()

    def show(self):
        plt.show()

    def kinship(self, K, nclusters=1, style=None):
        r"""Plot Kinship matrix.

        Parameters
        ----------
        K : array_like
            Kinship matrix.
        nclusters : int or str
            Number of blocks to be seen from the heatmap. It defaults to ``1``,
            which means that no ordering is performed. Pass 'auto' to
            automatically determine the number of clusters. Pass an integer to
            select the number of clusters.
        style : dict
            Keyword arguments forwarded to the
            :func:`matplotlib.axes.Axes.imshow` function.

        Examples
        --------
        .. plot::

            import numpy as np
            from matplotlib import pyplot as plt
            from limix.io.examples import numpy_kinship_file_example
            import limix

            K = np.load(numpy_kinship_file_example())
            p = limix.plot.get()
            p.kinship(K)
            p.show()
        """
        from .kinship import plot_kinship
        plot_kinship(K, nclusters, style, self.axes)

    def boxplot(self, df, style=None):
        r"""Box plot of variable values from different categories.

        Parameters
        ----------
        df : data_frame
            A data frame containing `value`, `variable`, and `category`
            columns.
        style : dict
            Keyword arguments forwarded to :func:`seaborn.boxplot`
            function.

        Examples
        --------
        .. plot::
            :include-source:

            import seaborn as sns
            import limix

            df = sns.load_dataset("exercise")
            df.rename(
                columns=dict(time='category', kind='variable', pulse='value'),
                inplace=True)

            p = limix.plot.get()
            p.boxplot(df)
            p.show()
        """
        from seaborn import boxplot

        if style is None:
            style = dict()

        boxplot(
            y='value',
            x='category',
            hue='variable',
            data=df,
            flierprops=dict(markersize=3.0),
            **style)

        self.axes.grid(
            True,
            which='major',
            axis='y',
            linewidth=0.75,
            linestyle='-',
            color='#EEEEEE',
            alpha=1.0)

    def curve(self, data, style=None):
        r"""Plot a curve and a confidence band around it.

        Parameters
        ----------
        data : list
            List of curves.
        style : dict
            Keyword arguments forwarded to :func:`matplotlib.axes.Axes.plot`
            function.

        Examples
        --------
        .. plot::

            from numpy.random import RandomState
            from numpy import sort
            import limix

            random = RandomState(0)

            x = random.randn(100)
            x = sort(x)
            y = x + random.randn(100) * 0.1
            ytop = y + 0.5
            ybottom = y - 0.5

            data = [dict(label='A', x=x, y=y, ybottom=ybottom,
                         ytop=ytop)]
            p = limix.plot.get()
            p.curve(data)
            p.show()
        """
        from .curve import plot_curve
        plot_curve(data, style, self.axes)

    def manhattan(self, df, alpha=None, null_style=None, alt_style=None):
        r"""Produce a manhattan plot.

        Parameters
        ----------
        df : :class:`pandas.DataFrame`
            A Pandas DataFrame containing columns pv for p-values, pos for
            base-pair positions, and chrom for chromossome names.
        alpha : float
            Threshold for significance. Defaults to ``0.01`` significance level
            (bonferroni-adjusted).
        null_style : dict
            Keyword arguments forwarded to the
            :func:`matplotlib.axes.Axes.plot` function when plotting the
            non-significant results.
        alt_style : dict
            Keyword arguments forwarded to the
            :func:`matplotlib.axes.Axes.plot` function when plotting the
            significant results.

        Examples
        --------
        .. plot::

            from numpy.random import RandomState
            from numpy import arange, ones, kron
            from pandas import DataFrame
            import limix

            random = RandomState(1)
            pv = random.rand(5000)
            pv[1200:1250] = random.rand(50)**4
            chrom  = kron(arange(1, 6), ones(1000))
            pos = kron(ones(5), arange(1, 1001))
            df = DataFrame(data=dict(pv=pv, chrom=chrom, pos=pos))
            p = limix.plot.get()
            p.manhattan(df)
            p.show()
        """
        from .manhattan import plot_manhattan
        plot_manhattan(df, alpha, null_style, alt_style, self.axes)

    def normal(self, x, bins=20, nstd=2, style=None):
        r"""Plot a fit of a normal distribution to the data in x.

        Parameters
        ----------
        x : array_like
            Values to be fitted.
        bins : int
            Number of histogram bins.
        nstd : float)
            Standard deviation multiplier.
        style : dict
            Keyword arguments forwarded to the :func:`matplotlib.axes.Axes.plot`
            function.

        Examples
        --------
        .. plot::

            import limix
            from numpy.random import RandomState

            random = RandomState(10)
            x = random.randn(100)
            p = limix.plot.get()
            p.normal(x, nstd=2)
            p.show()
        """
        from .normal import plot_normal
        plot_normal(x, bins, nstd, style, self.axes)

    def pca(self, X, style=None):
        r"""Plot the first two principal components of a design matrix.

        Parameters
        ----------
        X : array_like
            Design matrix.
        style : dict
            Keyword arguments forwarded to the :func:`matplotlib.pyplt.scatter`
            function.

        Examples
        --------
        .. plot::

            from numpy.random import RandomState
            import limix

            random = RandomState(0)
            X = random.randn(30, 10)

            p = limix.plot.get()
            p.pca(X)
            p.show()
        """
        from .pca import plot_pca
        plot_pca(X, style, self.axes)

    def power(self, df, style=None):
        r"""Plot number of hits across significance levels.

        Parameters
        ----------
        df : :class:`pandas.DataFrame`
            Data frame with `pv` and `label` columns.
        style : dict
            Keyword arguments forwarded to :func:`matplotlib.axes.Axes.plot`
            function.

        Examples
        --------
        .. plot::

            import limix
            from pandas import DataFrame
            from numpy.random import RandomState

            random = RandomState(1)
            nsnps = 10000

            pv0 = list(random.rand(nsnps))
            pv1 = list(0.7 * random.rand(nsnps))

            data = dict(pv=pv0 + pv1,
                        label=['label0'] * nsnps + ['label1'] * nsnps)
            df = DataFrame(data=data)
            p = limix.plot.get()
            p.power(df)
            p.show()
        """
        from .power import plot_power
        plot_power(df, style, self.axes)

    def qqplot(self, df, alpha=0.05, style=None):
        r"""Quantile-Quantile of observed p-values versus theoretical ones.

        Parameters
        ----------
        df : :class:`pandas.DataFrame`
            Data frame.
        alpha : float
            Significance level defining the band boundary.
            Defaults to ``0.05``.
        style : dict
            Keyword arguments forwarded to :func:`matplotlib.axes.Axes.plot`
            function.

        Returns
        -------
        :class:`matplotlib.axes.Axes`
            Axes.

        Examples
        --------
        .. plot::

            from pandas import DataFrame
            import limix
            from numpy.random import RandomState

            random = RandomState(1)

            pv0 = random.rand(10000)
            pv1 = random.rand(10000)

            pv = list(pv0) + list(pv1)
            label = ['label0'] * len(pv0) + ['label1'] * len(pv1)
            df = DataFrame(data=dict(pv=pv, label=label))
            p = limix.plot.get()
            p.qqplot(df)
        """
        from .qqplot import qqplot
        qqplot(df, alpha, style, self.axes)


def _set_rc_params():

    font = {
        'font.size': 10,
        'font.family': 'sans-serif',
        'font.sans-serif': ['Open Sans', 'Verdana', 'Arial', 'sans-serif'],
        'axes.labelsize': 14,
        'axes.titlesize': 16,
        'xtick.labelsize': 14,
        'ytick.labelsize': 14,
        'legend.fontsize': 14,
    }

    rcParams = {
        'figure.figsize': "5, 5",
        'axes.axisbelow': True,
        'lines.linewidth': 1.0,
        'patch.linewidth': 0.5,
        'lines.markersize': 2.0,
        'lines.markeredgewidth': 0.5,
        'xtick.major.width': 0.6,
        'ytick.major.width': 0.6,
        'xtick.minor.width': 0.3,
        'ytick.minor.width': 0.3,
        'xtick.major.pad': 2.5,
        'ytick.major.pad': 2.5,
        'xtick.direction': 'in',
        'ytick.direction': 'in',
        'savefig.transparent': True,
        'savefig.bbox': 'tight'
    }
    rcParams.update(font)
    mpl.rcParams.update(rcParams)
