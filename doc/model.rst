*****
Model
*****

We also export the underlying inference methods limix uses.  Its classes can be accessed
via limix module :mod:`limix.model` or via the package :mod:`glimix_core` itself. Both
ways should work identically.

.. currentmodule:: limix.model

Linear Mixed Models
===================

.. autoclass:: limix.model.lmm.LMM
    :members:

.. autoclass:: limix.model.lmm.Kron2Sum
    :members:

Generalised Linear Mixed Models
===============================

.. autoclass:: limix.model.glmm.GLMMNormal
    :members:

.. autoclass:: limix.model.glmm.GLMMExpFam
    :members:

Gaussian Process
================

.. autoclass:: limix.model.gp.GP
    :members:

Generalised Gaussian Process
============================

.. autoclass:: limix.model.ggp.ExpFamGP
    :members:

Covariance
==========

.. autoclass:: limix.model.cov.EyeCov
    :members:

.. autoclass:: limix.model.cov.FreeFormCov
    :members:

.. autoclass:: limix.model.cov.GivenCov
    :members:

.. autoclass:: limix.model.cov.LinearCov
    :members:

.. autoclass:: limix.model.cov.SumCov
    :members:

.. autoclass:: limix.model.cov.LRFreeFormCov
    :members:

.. autoclass:: limix.model.cov.Kron2SumCov
    :members:

Mean
====

.. autoclass:: limix.model.mean.OffsetMean
    :members:

.. autoclass:: limix.model.mean.LinearMean
    :members:

.. autoclass:: limix.model.mean.SumMean
    :members:

.. autoclass:: limix.model.mean.KronMean
    :members:
