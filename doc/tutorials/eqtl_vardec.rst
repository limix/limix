.. _ipython_directive:

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

Importing limix
---------------

.. ipython::

    In [1]: import limix

Downloading the data
--------------------

We are going to use a HDF5 file containg phenotype and genotypedd.

.. ipython::

    In [2]: url = "http://rest.s3for.me/limix/smith08.hdf5.bz2"

    In [3]: limix.download(url, verbose=False)

    @doctest
    In [4]: print(limix.filehash("smith08.hdf5.bz2"))
    cbad3ba229a117495453379000e4115cc70c2dac5776953bc15382e5ee9a71a6

    In [5]: limix.extract("smith08.hdf5.bz2", verbose=False)

    @doctest
    In [6]: limix.io.hdf5.see_hdf5("smith08.hdf5", verbose=False)
    /
      +--genotype
      |  +--col_header
      |  |  +--chrom [int64, (2956,)]
      |  |  +--pos [int64, (2956,)]
      |  |  +--pos_cum [int64, (2956,)]
      |  +--matrix [float64, (109, 2956)]
      |  +--row_header
      |     +--sample_ID [int64, (109,)]
      +--phenotype
         +--col_header
         |  +--environment [float64, (10986,)]
         |  +--gene_ID [|S100, (10986,)]
         |  +--gene_chrom [|S100, (10986,)]
         |  +--gene_end [int64, (10986,)]
         |  +--gene_start [int64, (10986,)]
         |  +--gene_strand [|S100, (10986,)]
         |  +--phenotype_ID [|S100, (10986,)]
         +--matrix [float64, (109, 10986)]
         +--row_header
            +--sample_ID [int64, (109,)]

    In [7]: data = limix.io.read_hdf5_limix("smith08.hdf5")

    @doctest
    In [8]: print(data['phenotype']['row_header'].head())
       sample_ID  i
    0          0  0
    1          1  1
    2          2  2
    3          3  3
    4          4  4

    @doctest
    In [9]: print(data['phenotype']['col_header'].head())
       environment  gene_ID gene_chrom  gene_end  gene_start gene_strand  \
    0          0.0  YOL161C         15     11548       11910           C
    1          0.0  YJR107W         10    628319      627333           W
    2          0.0  YPL270W         16     32803       30482           W
    3          0.0  YGR256W          7   1006108     1004630           W
    4          0.0  YDR518W          4   1480153     1478600           W

      phenotype_ID  i
    0    YOL161C:0  0
    1    YJR107W:0  1
    2    YPL270W:0  2
    3    YGR256W:0  3
    4    YDR518W:0  4

Selecting gene YBR115C under the glucose condition
--------------------------------------------------

.. ipython::

    In [10]: # Query for a specific phenotype

    In [11]: header = data['phenotype']['col_header']

    In [12]: query = "gene_ID=='YBR115C' and environment==0"

    In [13]: idx = header.query(query).i.values

    In [14]: # Select the phenotype itself

    In [15]: y = data['phenotype']['matrix'][:, idx].ravel()

    @savefig yeast_pheno01.png width=5in
    In [16]: limix.plot.plot_normal(y);

    In [17]: limix.plot.clf()

Genetic relatedness matrix
--------------------------

.. ipython:: python

    G = data['genotype']['matrix']
    K = limix.stats.linear_kinship(G, verbose=False)
    @savefig yeast_K01.png width=5in
    limix.plot.plot_kinship(K);

Univariate association test with linear mixed model
---------------------------------------------------

.. ipython:: python

    >>> result = limix.qtl.scan(G, y, 'normal', K, verbose=False)
    @doctest
    >>> print(result)
        Variants
              effsizes  effsizes_se       pvalues
    count  2956.000000  2956.000000  2.956000e+03
    mean      0.129739     0.589186  5.605584e-01
    std       0.550630     0.114092  2.778524e-01
    min      -1.267119     0.414053  2.583310e-20
    25%      -0.230129     0.518686  3.339200e-01
    50%       0.071479     0.563135  5.610395e-01
    75%       0.449852     0.611174  8.007013e-01
    max       4.198421     0.963061  9.996668e-01

    Covariate effect sizes for the null model
       offset
     0.012073

Cleaning up
-----------

.. ipython:: python

    import os
    from glob import glob
    for f in glob("smith08.hdf5*"):
        os.unlink(f)
