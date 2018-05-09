eQTL
^^^^

This tutorial illustrates the use of limix to analyse expression datasets.
For this illustration, we consider gene expression levels from a yeast genetics
study with freely available data . These data span 109 individuals with 2,956
marker SNPs and expression levels for 5,493 in glucose and ethanol growth media
respectively.
It is based on the `eQTL basics tutorial`_ of limix 1.0, which is now
deprecated.

.. _eQTL basics tutorial: https://github.com/limix/limix-tutorials/blob/master/eQTL/eQTL_basics.ipynb

Importing limix
---------------

.. nbplot::

    >>> import limix
    >>> import matplotlib.pyplot as plt

Downloading data
----------------

We are going to use a HDF5 file containg phenotype and genotyped data from
a remote repository.
Limix provides some handy utilities to perform common command line tasks,
like as downloading and extracting files.
However, feel free to use whatever method you prefer.

.. nbplot::

    >>> url = "http://rest.s3for.me/limix/smith08.hdf5.bz2"
    >>> limix.download(url, verbose=False)
    >>> print(limix.filehash("smith08.hdf5.bz2"))
    aecd5ebabd13ed2e38419c11d116e8d582077212efb37871a50c3a08fadb2ee1
    >>> limix.extract("smith08.hdf5.bz2", verbose=False)
    >>> print(limix.filehash("smith08.hdf5"))
    4648f596249ee2e3e60b9cd024d6f1af257079d39b2bff5192407a30de989266
    >>> limix.io.hdf5.see_hdf5("smith08.hdf5", verbose=False)
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
    >>> data = limix.io.read_hdf5_limix("smith08.hdf5")
    >>> print(data['phenotype']['row_header'].head())
       sample_ID  i
    0          0  0
    1          1  1
    2          2  2
    3          3  3
    4          4  4
    >>> print(data['phenotype']['col_header'].head())
       environment  gene_ID gene_chrom  gene_end  gene_start gene_strand  \
    0          0.0  YOL161C         15     11548       11910           C
    1          0.0  YJR107W         10    628319      627333           W
    2          0.0  YPL270W         16     32803       30482           W
    3          0.0  YGR256W          7   1006108     1004630           W
    4          0.0  YDR518W          4   1480153     1478600           W
    <BLANKLINE>
      phenotype_ID  i
    0    YOL161C:0  0
    1    YJR107W:0  1
    2    YPL270W:0  2
    3    YGR256W:0  3
    4    YDR518W:0  4

Selecting gene YBR115C under the glucose condition
--------------------------------------------------

Query for a specific phenotype, select the phenotype itself, and plot it.
The glucose condition is given by the environment ``0``.

.. nbplot::

    >>> header = data['phenotype']['col_header']
    >>> query = "gene_ID=='YBR115C' and environment==0"
    >>> idx = header.query(query).i.values
    >>> y = data['phenotype']['matrix'][:, idx].ravel()
    >>> limix.plot.normal(y)

Genetic relatedness matrix
--------------------------

The genetic relatedness will be determined by the inner-product of SNP
readings between individuals, and the result will be visualised via heatmap.

.. nbplot::

    >>> G = data['genotype']['matrix']
    >>> K = limix.stats.linear_kinship(G, verbose=False)
    >>> plt.clf()
    >>> limix.plot.kinship(K)

Univariate association test with linear mixed model
---------------------------------------------------

You have the option to either pass a raw array of samples-by-candidates for
the association scan or pass a tabular structure naming those candidates.
We recommend the second option as it will help maintain the association between
the results and the corresponding candidates.

The naming of those candidates is defined here by concatenating the chromossome
name and base-pair position.
However, it is often the case that SNP IDs are provided along with the
data, which can naturally be used for naming those candidates.

.. nbplot::

    >>> print(data['genotype']['col_header'].head())
    chrom   pos  pos_cum  i
    0      1   483      483  0
    1      1   484      484  1
    2      1  3220     3220  2
    3      1  3223     3223  3
    4      1  3232     3232  4
    >>> from pandas import DataFrame
    >>> chrom = data['genotype']['col_header']['chrom']
    >>> pos = data['genotype']['col_header']['pos']
    >>> candidate_ids = ["c{}_p{}".format(c, p) for c, p in zip(chrom, pos)]
    >>> G = DataFrame(G, columns=candidate_ids)
    >>> print(G.head())
    c1_p483  c1_p484  c1_p3220  c1_p3223  c1_p3232  c1_p3235  c1_p3244  \
    0      1.0      1.0       1.0       1.0       1.0       1.0       1.0
    1      1.0      0.0       1.0       1.0       1.0       1.0       1.0
    2      0.0      0.0       0.0       0.0       0.0       0.0       0.0
    3      0.0      0.0       1.0       1.0       1.0       1.0       1.0
    4      0.0      0.0       0.0       0.0       0.0       0.0       0.0
    <BLANKLINE>
    c1_p3247  c1_p3250  c1_p3274     ...       c16_p890898  c16_p890904  \
    0       1.0       1.0       1.0     ...               0.0          0.0
    1       1.0       1.0       1.0     ...               0.0          0.0
    2       0.0       0.0       0.0     ...               0.0          0.0
    3       1.0       1.0       1.0     ...               0.0          0.0
    4       0.0       0.0       0.0     ...               1.0          1.0
    <BLANKLINE>
    c16_p896709  c16_p897526  c16_p927500  c16_p927502  c16_p927506  \
    0          0.0          0.0          0.0          0.0          0.0
    1          0.0          0.0          1.0          1.0          1.0
    2          0.0          0.0          0.0          0.0          0.0
    3          0.0          0.0          0.0          0.0          0.0
    4          1.0          1.0          0.0          0.0          0.0
    <BLANKLINE>
    c16_p932310  c16_p932535  c16_p932538
    0          0.0          0.0          0.0
    1          1.0          1.0          1.0
    2          0.0          0.0          0.0
    3          0.0          1.0          1.0
    4          0.0          0.0          0.0
    <BLANKLINE>
    [5 rows x 2956 columns]

As you can see, we now have a pandas data frame ``G`` that keeps the candidate
identifications together with the actual allele read.
This data frame can be readily used to perform association scan.

.. nbplot::

    >>> qtl = limix.qtl.scan(G, y, 'normal', K, verbose=False)
    >>> print(qtl)
    Variants
          effsizes  effsizes_se       pvalues
    count  2956.000000  2956.000000  2.956000e+03
    mean      0.129739     0.589186  5.605584e-01
    std       0.550630     0.114092  2.778524e-01
    min      -1.267119     0.414053  2.583307e-20
    25%      -0.230129     0.518686  3.339200e-01
    50%       0.071479     0.563135  5.610395e-01
    75%       0.449852     0.611174  8.007013e-01
    max       4.198421     0.963061  9.996669e-01
    <BLANKLINE>
    Covariate effect sizes for the null model
    offset
    0.012073

Printing the result of an association scan will show a summary of the results.

Inspecting the p-values and effect-sizes are now easier because candidate
names are kept together with their corresponding statistics.

.. nbplot::

    >>> sorted_pvs = qtl.variant_pvalues.sort_values()
    >>> print(sorted_pvs.head())
    c2_p477206    2.583307e-20
    c2_p479161    1.250239e-13
    c2_p479164    1.250239e-13
    c2_p479166    1.250239e-13
    c2_p480009    9.086078e-13
    dtype: float64
    >>> print(qtl.variant_effsizes.loc[sorted_pvs.index].head())
    c2_p477206    4.198421
    c2_p479161    3.839388
    c2_p479164    3.839388
    c2_p479166    3.839388
    c2_p480009    3.857026
    dtype: float64

A Manhattan plot now automaticallt tags the significant associations using
their names.

.. nbplot::

    >>> pvs = qtl.variant_pvalues
    >>> pv = pvs.values
    >>> chrom = [i.split('_')[0][1:] for i, _ in pvs.iteritems()]
    >>> pos = [int(i.split('_')[1][1:]) for i, _ in pvs.iteritems()]
    >>> label = pvs.index.values
    >>> df = DataFrame(data=dict(pv=pv, chr=chrom, pos=pos, label=label))
    >>> plt.clf()
    >>> limix.plot.manhattan(df);
