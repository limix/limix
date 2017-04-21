r"""
***************************
Interaction Set Test (iSet)
***************************

- :func:`.fit_iSet`

.. rubric:: References

.. [CHRS17] Casale FP, Horta D, Rakitsch B, Stegle O (2017) Joint genetic analysis using variant sets reveals polygenic gene-context interactions. PLOS Genetics 13(4): e1006693.

Public interface
^^^^^^^^^^^^^^^^
"""

from .iset import fit_iSet

__all__ = ['fit_iSet']
