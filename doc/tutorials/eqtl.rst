eQTL
^^^^

This tutorial illustrates the use of limix to analyse expression datasets.
We consider gene expression levels from a yeast genetics
study with freely available data.
This data set span 109 individuals with 2,956 marker SNPs and expression
levels for 5,493 in glucose and ethanol growth media respectively.
It is based on the `eQTL basics tutorial`_ of limix 1.0, which is now
deprecated.

.. _eQTL basics tutorial: https://github.com/limix/limix-tutorials/blob/master/eQTL/eQTL_basics.ipynb

Importing limix
---------------

.. nbplot::

    >>> import limix

Downloading data
----------------

We are going to use a HDF5 file containg phenotype and genotyped data from
a remote repository.
Limix provides some handy utilities to perform common command line tasks,
like as downloading and extracting files.
However, feel free to use whatever method you prefer.

.. nbplot::

    >>> url = "http://rest.s3for.me/limix/smith08.hdf5.bz2"
    >>> limix.sh.download(url, verbose=False)
    >>> print(limix.sh.filehash("smith08.hdf5.bz2"))
    aecd5ebabd13ed2e38419c11d116e8d582077212efb37871a50c3a08fadb2ee1
    >>> _ = limix.sh.extract("smith08.hdf5.bz2", verbose=False)
    >>> limix.io.hdf5.see("smith08.hdf5")
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
    >>> data = limix.io.hdf5.read_limix("smith08.hdf5")
    >>> print(data['phenotype']['row_header'].head())
       sample_ID  i
    0          0  0
    1          1  1
    2          2  2
    3          3  3
    4          4  4
    >>> print(data['phenotype']['col_header'].head())
       environment  gene_ID gene_chrom  gene_end  gene_start gene_strand phenotype_ID  i
    0      0.00000  YOL161C         15     11548       11910           C    YOL161C:0  0
    1      0.00000  YJR107W         10    628319      627333           W    YJR107W:0  1
    2      0.00000  YPL270W         16     32803       30482           W    YPL270W:0  2
    3      0.00000  YGR256W          7   1006108     1004630           W    YGR256W:0  3
    4      0.00000  YDR518W          4   1480153     1478600           W    YDR518W:0  4

Selecting gene YBR115C under the glucose condition
--------------------------------------------------

Query for a specific phenotype, select the phenotype itself, and plot it.
The glucose condition is given by the environment ``0``.

.. nbplot::

    >>> header = data['phenotype']['col_header']
    >>> query = "gene_ID=='YBR115C' and environment==0"
    >>> idx = header.query(query).i.values
    >>> y = data['phenotype']['matrix'][:, idx].ravel()
    >>> limix.plot.normal(y) # doctest: +SKIP

Genetic relatedness matrix
--------------------------

The genetic relatedness will be determined by the inner-product of SNP
readings between individuals, and the result will be visualised via heatmap.

.. nbplot::

    >>> G = data['genotype']['matrix']
    >>> K = limix.stats.linear_kinship(G, verbose=False)
    >>> limix.plot.kinship(K) # doctest: +SKIP

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

    >>> from pandas import DataFrame
    >>> import numpy as np
    >>>
    >>> print(data['genotype']['col_header'].head())
    chrom   pos  pos_cum  i
    0      1   483      483  0
    1      1   484      484  1
    2      1  3220     3220  2
    3      1  3223     3223  3
    4      1  3232     3232  4
    >>> chrom = data['genotype']['col_header']['chrom']
    >>> pos = data['genotype']['col_header']['pos']
    >>> candidate_ids = ["c{}_p{}".format(c, p) for c, p in zip(chrom, pos)]
    >>> G = DataFrame(G, columns=candidate_ids)
    >>> print(G.head())
       c1_p483  c1_p484  c1_p3220  c1_p3223  c1_p3232  c1_p3235  c1_p3244  c1_p3247  \
    0  1.00000  1.00000   1.00000   1.00000   1.00000   1.00000   1.00000   1.00000
    1  1.00000  0.00000   1.00000   1.00000   1.00000   1.00000   1.00000   1.00000
    2  0.00000  0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000
    3  0.00000  0.00000   1.00000   1.00000   1.00000   1.00000   1.00000   1.00000
    4  0.00000  0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000
    <BLANKLINE>
       c1_p3250  c1_p3274  c1_p3280  c1_p3283  c1_p7292  c1_p7298  c1_p7358  c1_p7400  \
    0   1.00000   1.00000   1.00000   1.00000   1.00000   1.00000   1.00000   1.00000
    1   1.00000   1.00000   1.00000   1.00000   1.00000   1.00000   1.00000   1.00000
    2   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000
    3   1.00000   1.00000   1.00000   1.00000   1.00000   1.00000   1.00000   1.00000
    4   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000
    <BLANKLINE>
       c1_p7472  c1_p7478  c1_p7490  c1_p7532  c1_p7544  c1_p7574  c1_p7640  c1_p7652  \
    0   1.00000   1.00000   1.00000   1.00000   1.00000   1.00000   1.00000   1.00000
    1   1.00000   1.00000   1.00000   1.00000   1.00000   1.00000   1.00000   1.00000
    2   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000
    3   1.00000   1.00000   1.00000   1.00000   1.00000   1.00000   1.00000   1.00000
    4   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000
    <BLANKLINE>
       c1_p7712  c1_p10131  c1_p10134  c1_p10143  c1_p10146  c1_p10152  c1_p10236  \
    0   1.00000    1.00000    1.00000    1.00000    1.00000    1.00000    1.00000
    1   1.00000    1.00000    1.00000    1.00000    1.00000    1.00000    1.00000
    2   0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000
    3   1.00000    1.00000    1.00000    1.00000    1.00000    1.00000    1.00000
    4   0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000
    <BLANKLINE>
       c1_p10239  c1_p10284  c1_p10296  c1_p10302  c1_p10386  c1_p11582  c1_p11586  \
    0    1.00000    1.00000    1.00000    1.00000    1.00000    1.00000    1.00000
    1    1.00000    1.00000    1.00000    1.00000    1.00000    1.00000    1.00000
    2    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000
    3    1.00000    1.00000    1.00000    1.00000    1.00000    1.00000    1.00000
    4    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000
    <BLANKLINE>
       c1_p11588     ...       c16_p533282  c16_p535973  c16_p535979  c16_p542295  \
    0    1.00000     ...           1.00000      1.00000      1.00000      1.00000
    1    1.00000     ...           1.00000      1.00000      1.00000      1.00000
    2    0.00000     ...           1.00000      1.00000      1.00000      1.00000
    3    1.00000     ...           0.00000      0.00000      0.00000      0.00000
    4    0.00000     ...           0.00000      0.00000      0.00000      0.00000
    <BLANKLINE>
       c16_p542307  c16_p547618  c16_p555416  c16_p590622  c16_p600658  c16_p600664  \
    0      1.00000      1.00000      1.00000      0.00000      0.00000      0.00000
    1      1.00000      1.00000      1.00000      1.00000      1.00000      1.00000
    2      1.00000      1.00000      1.00000      1.00000      1.00000      1.00000
    3      0.00000      0.00000      0.00000      0.00000      0.00000      0.00000
    4      0.00000      0.00000      0.00000      0.00000      0.00000      0.00000
    <BLANKLINE>
       c16_p604010  c16_p618575  c16_p618581  c16_p620596  c16_p695782  c16_p700280  \
    0      0.00000      0.00000      0.00000      0.00000      0.00000      0.00000
    1      1.00000      1.00000      1.00000      1.00000      1.00000      1.00000
    2      1.00000      1.00000      1.00000      1.00000      1.00000      1.00000
    3      0.00000      0.00000      0.00000      0.00000      0.00000      1.00000
    4      0.00000      0.00000      0.00000      0.00000      0.00000      0.00000
    <BLANKLINE>
       c16_p704388  c16_p711614  c16_p718892  c16_p718893  c16_p744530  c16_p744590  \
    0      0.00000      0.00000      0.00000      0.00000      0.00000      0.00000
    1      1.00000      1.00000      1.00000      1.00000      1.00000      1.00000
    2      1.00000      1.00000      1.00000      1.00000      1.00000      1.00000
    3      1.00000      1.00000      1.00000      1.00000      0.00000      0.00000
    4      0.00000      0.00000      1.00000      1.00000      1.00000      1.00000
    <BLANKLINE>
       c16_p744599  c16_p748158  c16_p787283  c16_p819247  c16_p819249  c16_p819251  \
    0      0.00000      0.00000      0.00000      0.00000      0.00000      0.00000
    1      1.00000      1.00000      1.00000      0.00000      0.00000      0.00000
    2      1.00000      1.00000      1.00000      1.00000      1.00000      1.00000
    3      0.00000      0.00000      0.00000      1.00000      1.00000      1.00000
    4      1.00000      1.00000      1.00000      1.00000      1.00000      1.00000
    <BLANKLINE>
       c16_p825431  c16_p890898  c16_p890904  c16_p896709  c16_p897526  c16_p927500  \
    0      0.00000      0.00000      0.00000      0.00000      0.00000      0.00000
    1      0.00000      0.00000      0.00000      0.00000      0.00000      1.00000
    2      0.00000      0.00000      0.00000      0.00000      0.00000      0.00000
    3      1.00000      0.00000      0.00000      0.00000      0.00000      0.00000
    4      1.00000      1.00000      1.00000      1.00000      1.00000      0.00000
    <BLANKLINE>
       c16_p927502  c16_p927506  c16_p932310  c16_p932535  c16_p932538
    0      0.00000      0.00000      0.00000      0.00000      0.00000
    1      1.00000      1.00000      1.00000      1.00000      1.00000
    2      0.00000      0.00000      0.00000      0.00000      0.00000
    3      0.00000      0.00000      0.00000      1.00000      1.00000
    4      0.00000      0.00000      0.00000      0.00000      0.00000
    <BLANKLINE>
    [5 rows x 2956 columns]

As you can see, we now have a pandas data frame ``G`` that keeps the candidate
identifications together with the actual allele read.
This data frame can be readily used to perform association scan.

.. nbplot::

    >>> qtl = limix.qtl.scan(G, y, 'normal', K, verbose=False)
    >>> print(qtl) # doctest: +FLOAT_CMP
    Variants
    --------
            effsizes  effsizes_se    pvalues
    count 2956.00000   2956.00000 2956.00000
    mean     0.12974      0.58919    0.56056
    std      0.55063      0.11409    0.27785
    min     -1.26712      0.41405    0.00000
    25%     -0.23013      0.51869    0.33392
    50%      0.07148      0.56313    0.56104
    75%      0.44985      0.61117    0.80070
    max      4.19842      0.96306    0.99967
    <BLANKLINE>
    Covariate effect sizes for H0
    -----------------------------
     0
        0.01207

Inspecting the p-values and effect-sizes are now easier because candidate
names are kept together with their corresponding statistics.

.. nbplot::

    >>> pv = qtl.variant_pvalues.sort_values()
    >>> print(np.log(pv.head())) # doctest: +FLOAT_CMP
    c2_p477206   -45.10263
    c2_p479161   -29.71027
    c2_p479164   -29.71027
    c2_p479166   -29.71027
    c2_p480009   -27.72686
    dtype: float64
    >>> print(qtl.variant_effsizes.loc[pv.index].head()) # doctest: +FLOAT_CMP
    c2_p477206    4.19842
    c2_p479161    3.83939
    c2_p479164    3.83939
    c2_p479166    3.83939
    c2_p480009    3.85703
    dtype: float64

A Manhattan plot can help understand the result.

.. nbplot::

    >>> pv = qtl.variant_pvalues
    >>> chrom = [i.split('_')[0][1:] for i, _ in pv.iteritems()]
    >>> pos = [int(i.split('_')[1][1:]) for i, _ in pv.iteritems()]
    >>> df = DataFrame(data=dict(pv=pv, chr=chrom, pos=pos))
    >>> limix.plot.manhattan(df) # doctest: +SKIP

We then remove the temporary files.

.. nbplot::

    >>> limix.sh.remove("smith08.hdf5.bz2")
    >>> limix.sh.remove("smith08.hdf5")
