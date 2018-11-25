r"""
Scanning multiple variants
^^^^^^^^^^^^^^^^^^^^^^^^^^

For multi-variant scans with StructLMM see :ref:`gwas`.

Interpretation tools
^^^^^^^^^^^^^^^^^^^^

- :class:`.OptimalRho`
- :class:`.PredictGenEffect`
- :class:`.BF`

"""

from ._optimal_rho import OptimalRho
from ._predict_genetic_effect import PredictGenEffect
from ._bf import BF

__all__ = ["OptimalRho", "PredictGenEffect", "BF"]
