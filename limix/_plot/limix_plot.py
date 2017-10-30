import matplotlib as mpl
from matplotlib import pyplot as plt

_figure_id = 98748


class LimixPlot(object):
    def __init__(self):
        self._figure = None
        self._axes = None
        self._figure_initialised = False

    def _initialise_figure(self):
        plt.close(_figure_id)

        _matplotlib_setup()
        figure = plt.figure(_figure_id)
        axes = figure.add_subplot(111)

        self._figure = figure
        self._axes = axes
        self._figure_initialised = True

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

    def savefig(self, *args, **kwargs):
        self._figure.savefig(*args, **kwargs)

    def show(self):
        plt.show()

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

        Example
        -------
        .. plot::
            :include-source:

            import limix

            df = limix.load_dataset("boxplot")

            limix.plot.boxplot(df)
            limix.plot.show()
        """
        from .boxplot import plot_boxplot
        self._initialise_figure()
        plot_boxplot(df, style, self.axes)

    def curve(self, data, style=None):
        r"""Plot a curve and a confidence band around it.

        Parameters
        ----------
        data : list
        List of curves.
        style : dict
        Keyword arguments forwarded to :func:`matplotlib.axes.Axes.plot`
        function.

        Example
        -------
        .. plot::
            :include-source:

            import limix
            from numpy.random import RandomState
            from numpy import sort

            random = RandomState(0)
            x = random.randn(100)
            x = sort(x)
            y = x + random.randn(100) * 0.1
            ytop = y + 0.5
            ybottom = y - 0.5
            data = [dict(label='A', x=x, y=y, ybottom=ybottom,
            ytop=ytop)]

            limix.plot.curve(data)
            limix.plot.show()
        """
        from .curve import plot_curve
        self._initialise_figure()
        plot_curve(data, style, self.axes)

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

        Example
        -------
        .. plot::
            :include-source:

            import limix

            K = limix.load_dataset("kinship", verbose=False)

            limix.plot.kinship(K)
            limix.plot.show()
        """
        from .kinship import plot_kinship
        self._initialise_figure()
        plot_kinship(K, nclusters, style, self.axes)

    def manhattan(self, df, alpha=None, null_style=None, alt_style=None):
        r"""Produce a manhattan plot.

        Parameters
        ----------
        df : :class:`pandas.DataFrame`
        A Pandas DataFrame containing columns `pv` for p-values, `pos` for
        base-pair positions, and `chrom` for chromossome names.
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

        Example
        -------
        .. plot::
            :include-source:

            import limix
            from numpy.random import RandomState
            from numpy import arange, ones, kron
            from pandas import DataFrame

            random = RandomState(1)
            pv = random.rand(5000)
            pv[1200:1250] = random.rand(50)**4
            chrom  = kron(arange(1, 6), ones(1000))
            pos = kron(ones(5), arange(1, 1001))
            df = DataFrame(data=dict(pv=pv, chrom=chrom, pos=pos))

            limix.plot.manhattan(df)
            limix.plot.show()
        """
        from .manhattan import plot_manhattan
        self._initialise_figure()
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

        Example
        -------
        .. plot::
            :include-source:

            import limix
            from numpy.random import RandomState

            random = RandomState(10)
            x = random.randn(100)
            limix.plot.normal(x, nstd=2)
            limix.plot.show()
        """
        from .normal import plot_normal
        self._initialise_figure()
        plot_normal(x, bins, nstd, style, self.axes)

    def pca(self, X, style=None):
        r"""Plot the first two principal components of a design matrix.

        Parameters
        ----------
        X : array_like
        Design matrix.
        style : dict
        Keyword arguments forwarded to the :func:`matplotlib.pyplt.plot`
        function.

        Example
        -------
        .. plot::
            :include-source:

            import limix
            from numpy.random import RandomState

            random = RandomState(0)
            X = random.randn(30, 10)

            limix.plot.pca(X)
            limix.plot.show()
        """
        from .pca import plot_pca
        self._initialise_figure()
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

        Example
        -------
        .. plot::
            :include-source:

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

            limix.plot.power(df)
            limix.plot.show()
        """
        from .power import plot_power
        self._initialise_figure()
        plot_power(df, style, self.axes)

    def qqplot(self, df, alpha=0.05, cutoff=0.1, style=None):
        r"""Quantile-Quantile of observed p-values versus theoretical ones.

        Parameters
        ----------
        df : :class:`pandas.DataFrame`
        Data frame.
        alpha : float
        Significance level defining the band boundary.
        Defaults to ``0.05``.
        cutoff : float
        P-values higher than `cutoff` will not be plotted.
        style : dict
        Keyword arguments forwarded to :func:`matplotlib.axes.Axes.plot`
        function.

        Example
        -------
        .. plot::
            :include-source:

            import limix
            from pandas import DataFrame
            from numpy.random import RandomState

            random = RandomState(1)

            pv0 = random.rand(10000)
            pv1 = random.rand(10000)

            pv = list(pv0) + list(pv1)
            label = ['label0'] * len(pv0) + ['label1'] * len(pv1)
            df = DataFrame(data=dict(pv=pv, label=label))

            limix.plot.qqplot(df)
            limix.plot.show()
        """
        from .qqplot import qqplot
        self._initialise_figure()
        qqplot(df, alpha, cutoff, style, self.axes)

    def see_image(self, filepath, verbose=True):
        r"""Plot a image represented in a file.

        It uses :func:`matplotlib.pyplot.imread` in order to read the provided
        file and plot.

        Parameters
        ----------
        filepath : str
        File path to the image.
        """
        from .image import see_image
        see_image(filepath, verbose)


def _matplotlib_setup():

    font = {
        'font.size': 10,
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
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
