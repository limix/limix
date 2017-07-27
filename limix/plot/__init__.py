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
from .power import plot_power, plot_power_known
from .qqplot import plot_qqplot


def show(*args, **kwargs):
    r"""Show a plot."""
    from matplotlib import pyplot as plt
    plt.tight_layout()
    plt.show(*args, **kwargs)


def figure(*args, **kwargs):
    r"""Create a figure."""
    from matplotlib import pyplot as plt
    return plt.figure(*args, **kwargs)

def savefig(*args, **kwargs):
    r"""Save a figure."""
    from matplotlib import pyplot as plt
    return plt.savefig(*args, **kwargs)


__all__ = [
    'plot_manhattan', 'plot_qqplot', 'plot_normal', 'plot_kinship',
    'see_image', 'plot_power', 'plot_power_known', 'plot_curve',
    'ConsensusCurve', 'show', 'figure', 'plot_qqnormal', 'savefig'
]
