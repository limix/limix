# TODO: document it
def set_paper_style():
    import matplotlib as mpl
    from seaborn import set_context, despine

    fonts = [
        'Helvetica', 'Arial', 'Baskerville', 'Caslon', 'Garamond',
        'DejaVu Sans', 'Bitstream Vera Sans', 'Computer Modern Sans Serif',
        'Lucida Grande', 'Verdana', 'Geneva', 'Lucid', 'Avant Garde',
        'sans-serif'
    ]

    mpl.pyplot.gca().set_axisbelow(True)
    despine(left=False, right=False, top=False, bottom=False)

    rcParams = {
        'font.size': 10,
        'font.sans-serif': fonts,
        'axes.labelsize': 12,
        'axes.titlesize': 8,
        'axes.axisbelow': True,
        'xtick.labelsize': 12,
        'ytick.labelsize': 12,
        'legend.fontsize': 8.0,
        'lines.linewidth': 1.0,
        'patch.linewidth': 0.2,
        'lines.markersize': 5.6,
        'lines.markeredgewidth': 0.0,
        'xtick.major.width': 0.6,
        'ytick.major.width': 0.6,
        'xtick.minor.width': 0.3,
        'ytick.minor.width': 0.3,
        'xtick.major.pad': 5.6,
        'ytick.major.pad': 5.6
    }
    mpl.rcParams.update(rcParams)
