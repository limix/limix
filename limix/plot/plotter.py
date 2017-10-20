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

    def curve(self):
        pass

    def kinship(self):
        pass

    def manhattan(self):
        pass

    def normal(self):
        pass

    def qqnormal(self):
        pass

    def power(self):
        pass

    def power_known(self):
        pass

    def qqplot(self):
        pass

    def boxplot(self, df, palette):
        from seaborn import boxplot

        boxplot(
            y='value',
            x='category',
            hue='variable',
            data=df,
            flierprops=dict(markersize=3.0),
            palette=palette)

        self.axes.grid(
            True,
            which='major',
            axis='y',
            linewidth=0.75,
            linestyle='-',
            color='#EEEEEE',
            alpha=1.0)

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

            import pandas as pd
            from limix import get_plot
            from numpy.random import RandomState

            random = RandomState(1)

            pv0 = random.rand(10000)
            pv1 = random.rand(10000)

            data = dict(pv=list(pv0) + list(pv1),
                        label=['label0'] * len(pv0) + ['label1'] * len(pv1))

            p = limix.get_plot()
            p.qqplot(pd.DataFrame(data=data))
            p.show()
        """
        from .qqplot import qqplot
        qqplot(df, self.axes, alpha, style)


def _set_rc_params():

    font = {
        'font.size': 10,
        'font.family': 'sans-serif',
        'font.sans-serif': ['Open Sans', 'Verdana', 'Arial', 'sans-serif'],
        'axes.labelsize': 10,
        'axes.titlesize': 14,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'legend.fontsize': 12,
    }

    rcParams = {
        'figure.figsize': "5, 5",
        'axes.axisbelow': True,
        'lines.linewidth': 1.0,
        'patch.linewidth': 0.2,
        'lines.markersize': 2.0,
        'lines.markeredgewidth': 0.0,
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
