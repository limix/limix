from __future__ import division

from numpy import argsort, asarray


def plot_curve(data, style=None, ax=None):
    r"""Plot a curve with a band.

    Parameters
    ----------
    data : dict
        TODO.
    style : dict
        Keyword arguments forwarded to :func:`matplotlib.axes.Axes.plot`
        function.
    ax : :class:`matplotlib.axes.Axes`:
        The target handle for this figure. If None, the current axes is set.

    Returns
    -------
    :class:`matplotlib.axes.Axes`
        Axes object.

    Examples
    --------
    .. plot::

        import matplotlib.pyplot as plt
        import numpy as np
        import pandas as pd
        from limix.plot import plot_curve

        random = np.random.RandomState(0)

        x = random.randn(100)
        x = np.sort(x)
        y = x + random.randn(100) * 0.1
        ytop = y + 0.5
        ybottom = y - 0.5

        data = [dict(label='methodA', x=x, y=y, ybottom=ybottom, ytop=ytop)]
        plot_curve(data)

        plt.show()
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
                alpha=0.1,
                interpolate=True)

            ax.fill_between(
                x,
                ybottom,
                ytop,
                lw=0.05,
                edgecolor=color,
                facecolor='None',
                interpolate=True)

    ax.legend()

    return ax


if __name__ == '__main__':
    import numpy as np
    import pandas as pd
    random = np.random.RandomState(0)

    x = random.randn(100)
    x = np.sort(x)
    y = x + random.randn(100) * 0.1
    ytop = y + 0.5
    ybottom = y - 0.5

    data = [dict(label='methodA', x=x, y=y, ybottom=ybottom, ytop=ytop)]
    plot_curve(data)
    import matplotlib.pyplot as plt
    plt.savefig('/mnt/hd0/Scratch/fig.png')
