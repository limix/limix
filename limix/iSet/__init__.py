r"""
***************************
Interaction Set Test (iSet)
***************************

Joint genetic models for multiple traits have helped to enhance association
analyses.
Most existing multi-trait models have been designed to increase power for
detecting associations, whereas the analysis of interactions has received
considerably less attention.
Here, we implement :func:`.fit_iSet` [CHRS17]_, a method based on linear
mixed models to test for interactions between sets of variants and
environmental states or other contexts.
Our model generalizes previous
interaction tests and in particular provides a test for local differences in
the genetic architecture between contexts.
We first use simulations to validate
iSet before applying the model to the analysis of genotype-environment
interactions in an eQTL study.
Our model retrieves a larger number of
interactions than alternative methods and reveals that up to 20%% of cases show
context-specific configurations of causal variants.
Finally, we apply iSet to
test for sub-group specific genetic effects in human lipid levels in a large
human cohort, where we identify a gene-sex interaction for C-reactive protein
that is missed by alternative methods.

.. rubric:: References

.. [CHRS17] Casale FP, Horta D, Rakitsch B, Stegle O (2017) Joint genetic analysis using variant sets reveals polygenic gene-context interactions. PLOS Genetics 13(4): e1006693.

Public interface
^^^^^^^^^^^^^^^^
"""

from .iset import fit_iSet

__all__ = ['fit_iSet']
