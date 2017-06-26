r"""
******************
Plotting utilities
******************

Manhattan plot
^^^^^^^^^^^^^^

.. autofunction:: limix.plot.plot_manhattan

QQ plot
^^^^^^^

.. autofunction:: limix.plot.plot_qqplot

Power plots
^^^^^^^^^^^

.. autofunction:: limix.plot.plot_power_curve

Kinship plot
^^^^^^^^^^^^

.. autofunction:: limix.plot.plot_kinship

Normal distribution
^^^^^^^^^^^^^^^^^^^

.. autofunction:: limix.plot.plot_normal

File plot
^^^^^^^^^

.. autofunction:: limix.plot.see_image

"""

from .image import see_image
from .kinship import plot_kinship
from .manhattan import plot_manhattan
from .normal import plot_normal
from .power import plot_power_curve
from .qqplot import plot_qqplot

__all__ = [
    'plot_manhattan', 'plot_qqplot', 'plot_normal', 'plot_kinship',
    'see_image', 'plot_power_curve'
]
