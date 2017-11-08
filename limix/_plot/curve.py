from __future__ import division


def plot_curve(data, style=None, ax=None):
    r"""Plot a curve with a band.

    Parameters
    ----------
    data : list
        List of curves.
    style : dict
        Keyword arguments forwarded to :func:`matplotlib.axes.Axes.plot`
        function.
    ax : :class:`matplotlib.axes.Axes`:
        The target handle for this figure. If None, the current axes is set.

    Returns
    -------
    :class:`matplotlib.axes.Axes`
        Axes object.
    """

    import matplotlib.pyplot as plt

    ax = plt.gca() if ax is None else ax

    labels = [d['label'] for d in data]
    if style is None:
        style = {label: dict() for label in labels}

    for d in data:

        label = d['label']
        x = d['x']
        y = d['y']

        lines = ax.plot(x, y, label=label, **style.get(label))

        if 'ybottom' in d:
            if 'color' in style[label]:
                color = style[label]['color']
            else:
                color = lines[0].get_color()

            ybottom = d['ybottom']
            ytop = d['ytop']

            ax.fill_between(
                x,
                ybottom,
                ytop,
                lw=0,
                edgecolor='None',
                facecolor=color,
                alpha=0.25,
                interpolate=True)

    ax.legend(loc='best')
    ax.grid(True, which='major', axis='both', alpha=1.0)

    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')

    return ax
