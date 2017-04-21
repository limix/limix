r"""
****************************
Multi Trait Set Test (mtSet)
****************************

Set tests are a powerful approach for genome-wide association testing between
groups of genetic variants and quantitative traits.
Here we implement :class:`.MTSet`, a mixed-model approach that enables joint
analysis across multiple correlated traits while accounting for population
structure and relatedness.
:class:`.MTSet` effectively combines the benefits of set tests with multi-trait
modeling and is computationally efficient, enabling genetic analysis of large
cohorts (up to 500,000 individuals) and multiple traits.
A detailed description of the method can be found at [CRS15]_.

.. rubric:: References

.. [CRS15] Casale FP, Rakitsch B, Lippert C, Stegle O (2015) Efficient set tests for the genetic analysis of correlated traits. Nature Methods, Vol. 12, No. 8. (15 June 2015), pp. 755-758.

Public interface
^^^^^^^^^^^^^^^^
"""

from .mtset import MTSet

__all__ = ['MTSet']
