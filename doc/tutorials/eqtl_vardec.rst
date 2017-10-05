eQTL and variance decomposition
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


This tutorial illustrates the use of LIMIX to anlayse expression datasets. For
this illustration, we consider gene expression levels from a yeast genetics
study with freely available data . These data span 109 individuals with 2,956
marker SNPs and expression levels for 5,493 in glucose and ethanol growth media
respectively. We start out by discussing how to do QTL mapping, implement models
that consider multi loci and introduce the application of variance component
models for single quantitative traits. Subsequently, these analysis are extended
to the corresponding multi-trait models.

# .. doctest::

#     >>> import limix

# .. doctest::

#     >>> limix.download("http://rest.s3for.me/limix/smith08.hdf5.bz2")
