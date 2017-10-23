# TODO: document it
def set_paper_style():
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    ax = plt.gca()
    ax.set_axisbelow(True)
    ax.spines['right'].set_visible(True)
    ax.spines['top'].set_visible(True)

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
