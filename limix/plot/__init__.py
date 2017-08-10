r"""
******************
Plotting utilities
******************

Curve
^^^^^

.. autofunction:: limix.plot.plot_curve
.. autoclass:: limix.plot.ConsensusCurve

Manhattan plot
^^^^^^^^^^^^^^

.. autofunction:: limix.plot.plot_manhattan

QQ plot
^^^^^^^

.. autofunction:: limix.plot.plot_qqplot

Power plots
^^^^^^^^^^^

.. autofunction:: limix.plot.plot_power

Principal components plot
^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: limix.plot.plot_pca

Kinship plot
^^^^^^^^^^^^

.. autofunction:: limix.plot.plot_kinship

Normal distribution
^^^^^^^^^^^^^^^^^^^

.. autofunction:: limix.plot.plot_normal
.. autofunction:: limix.plot.plot_qqnormal

File plot
^^^^^^^^^

.. autofunction:: limix.plot.see_image

"""

from .consensus import ConsensusCurve
from .curve import plot_curve
from .image import see_image
from .kinship import plot_kinship
from .manhattan import plot_manhattan
from .normal import plot_normal, plot_qqnormal
from .pca import plot_pca
from .power import plot_power, plot_power_known
from .qqplot import plot_qqplot
from .style import set_paper_style


def show(*args, **kwargs):
    r"""Show a plot."""
    from matplotlib import pyplot as plt
    plt.show(*args, **kwargs)


def figure(*args, **kwargs):
    r"""Create a figure."""
    from matplotlib import pyplot as plt
    return plt.figure(*args, **kwargs)


def savefig(*args, **kwargs):
    r"""Save a figure."""
    from matplotlib import pyplot as plt
    return plt.savefig(*args, **kwargs)


def clf(*args, **kwargs):
    r"""Clear current figure."""
    from matplotlib import pyplot as plt
    return plt.clf(*args, **kwargs)


def set_xlabel(name):
    from matplotlib import pyplot as plt
    plt.gca().set_xlabel(name)


def set_ylabel(name):
    from matplotlib import pyplot as plt
    plt.gca().set_ylabel(name)


def set_xlim(left, right):
    from matplotlib import pyplot as plt
    plt.gca().set_xlim(left, right)


def set_ylim(left, right):
    from matplotlib import pyplot as plt
    plt.gca().set_ylim(left, right)


def legend(*args, **kwargs):
    from matplotlib import pyplot as plt
    plt.gca().legend(*args, **kwargs)


def tight_layout(*args, **kwargs):
    from matplotlib import pyplot as plt
    plt.tight_layout(*args, **kwargs)


def spine(left=True, right=True, top=True, bottom=True):
    from seaborn import despine
    spine(left=not Left, right=not right, top=not top, bottom=not bottom)


__all__ = [
    'plot_manhattan', 'plot_qqplot', 'plot_normal', 'plot_kinship',
    'see_image', 'plot_power', 'plot_power_known', 'plot_curve',
    'ConsensusCurve', 'show', 'figure', 'plot_qqnormal', 'savefig', 'plot_pca',
    'clf', 'set_paper_style', 'set_xlabel', 'set_ylabel', 'legend', 'set_xlim',
    'set_ylim', 'tight_layout'
]
