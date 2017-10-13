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
