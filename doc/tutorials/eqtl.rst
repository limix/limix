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

.. plot::
    :context:

    >>> from __future__ import unicode_literals
    >>> import limix

Downloading data
----------------

We are going to use a HDF5 file containg phenotype and genotyped data from
a remote repository.
Limix provides some handy utilities to perform common command line tasks,
like as downloading and extracting files.
However, feel free to use whatever method you prefer.

.. plot::
    :context:

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
    >>> Y = data['phenotype']
    >>> G = data['genotype']
    >>> print(Y)  # doctest: +ALLOW_UNICODE
    <xarray.DataArray 'phenotype' (sample: 109, outcome: 10986)>
    array([[-0.037339, -0.078165,  0.042936, ...,  0.095596, -0.132385, -0.274954],
           [-0.301376,  0.066055,  0.338624, ..., -0.142661, -0.238349,  0.732752],
           [ 0.002661,  0.121835, -0.137064, ..., -0.144404,  0.257615,  0.015046],
           ...,
           [-0.287339,  0.351835,  0.072936, ...,  0.097339, -0.038349,  0.162752],
           [-0.577339,  0.011835, -0.007064, ...,  0.135596,  0.107615,  0.245046],
           [-0.277339,  0.061835,  0.132936, ...,  0.015596, -0.142385, -0.124954]])
    Coordinates:
      * sample        (sample) int64 0 1 2 3 4 5 6 7 ... 102 103 104 105 106 107 108
        environment   (outcome) float64 0.0 0.0 0.0 0.0 0.0 ... 1.0 1.0 1.0 1.0 1.0
        gene_ID       (outcome) object 'YOL161C' 'YJR107W' ... 'YLR118C' 'YBR242W'
        gene_chrom    (outcome) object '15' '10' '16' '7' '4' ... '3' '10' '12' '2'
        gene_end      (outcome) int64 11548 628319 32803 ... 315049 384726 705381
        gene_start    (outcome) int64 11910 627333 30482 ... 315552 385409 704665
        gene_strand   (outcome) object 'C' 'W' 'W' 'W' 'W' ... 'W' 'W' 'C' 'C' 'W'
        phenotype_ID  (outcome) object 'YOL161C:0' 'YJR107W:0' ... 'YBR242W:1'
    Dimensions without coordinates: outcome

Selecting gene YBR115C under the glucose condition
--------------------------------------------------

Query for a specific phenotype, select the phenotype itself, and plot it.
The glucose condition is given by the environment ``0``.

header = data['phenotype']['col_header']
query = "gene_ID=='YBR115C' and environment==0"
idx = header.query(query).i.values
y = data['phenotype']['matrix'][:, idx].ravel()

.. plot::
    :context:

    >>> y = Y[:, (Y.gene_ID == "YBR115C") & (Y.environment==0)]
    >>> y = y.stack(z=('sample', 'outcome')).reset_index('z')
    >>> y = y.rename(z="sample")
    >>> _ = limix.plot.normal(y)
    >>> limix.plot.show()

Genetic relatedness matrix
--------------------------

The genetic relatedness will be determined by the inner-product of SNP
readings between individuals, and the result will be visualised via heatmap.

.. plot::
    :context:

    >>> K = limix.stats.linear_kinship(G.values, verbose=False)
    >>> _ = limix.plot.kinship(K)
    >>> limix.plot.show()

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

.. plot::
    :context:

    >>> from pandas import DataFrame
    >>> import numpy as np
    >>>
    >>> print(G)
    <xarray.DataArray 'genotype' (sample: 109, candidate: 2956)>
    array([[1., 1., 1., ..., 0., 0., 0.],
           [1., 0., 1., ..., 1., 1., 1.],
           [0., 0., 0., ..., 0., 0., 0.],
           ...,
           [0., 0., 0., ..., 0., 1., 1.],
           [0., 0., 0., ..., 1., 1., 1.],
           [1., 1., 1., ..., 1., 1., 1.]])
    Coordinates:
      * sample   (sample) int64 0 1 2 3 4 5 6 7 ... 101 102 103 104 105 106 107 108
        chrom    (candidate) int64 1 1 1 1 1 1 1 1 1 ... 16 16 16 16 16 16 16 16 16
        pos      (candidate) int64 483 484 3220 3223 ... 927506 932310 932535 932538
        pos_cum  (candidate) int64 483 484 3220 3223 ... 12055570 12055795 12055798
    Dimensions without coordinates: candidate

As you can see, we now have a pandas data frame ``G`` that keeps the candidate
identifications together with the actual allele read.
This data frame can be readily used to perform association scan.

.. plot::
    :context:

    >>> print(y)
    <xarray.DataArray 'phenotype' (sample: 109)>
    array([ 3.504479,  1.914585,  3.434479, -2.075521,  1.654585,  3.304479,
            2.044585, -4.125415,  2.024585,  0.732574, -0.80732 ,  3.464479,
           -2.385521,  3.644479, -4.785415, -2.895521,  0.732574, -3.155521,
           -0.80732 , -0.80732 ,  0.732574, -2.695521, -2.835521, -4.635415,
            0.732574,  1.804585, -0.80732 ,  1.964585,  2.304585,  2.484585,
            2.424585,  2.534585,  3.254479, -0.80732 , -4.555415, -1.815521,
            1.934585, -2.065521,  1.754585,  2.014585,  0.732574, -2.835521,
           -2.715521, -3.115521,  1.854585,  3.544479, -0.80732 ,  2.594585,
            3.574479, -3.175521,  0.732574, -2.395521,  1.824585,  2.134585,
           -0.80732 , -2.775521, -2.255521, -0.80732 ,  3.544479,  1.894585,
            3.364479, -2.775521, -0.80732 ,  2.034585, -4.695415, -0.80732 ,
            1.864585,  2.174585, -3.815521,  1.674585, -2.725521, -2.685521,
           -1.345521, -2.405521, -4.035415, -0.80732 , -0.80732 , -2.255521,
           -2.765521,  3.314479, -0.80732 ,  3.594479, -2.815521,  3.954479,
            1.794585,  1.904585,  2.064585, -0.80732 ,  3.864479,  3.604479,
           -2.505521, -0.80732 ,  1.804585, -2.345521, -0.80732 , -3.135521,
            3.704479,  3.714479, -4.565415,  0.732574, -0.80732 ,  0.732574,
            2.244585, -2.385521,  3.304479, -2.895521, -2.475521, -2.625521,
            3.314479])
    Coordinates:
        environment   (sample) float64 0.0 0.0 0.0 0.0 0.0 ... 0.0 0.0 0.0 0.0 0.0
        gene_ID       (sample) object 'YBR115C' 'YBR115C' ... 'YBR115C' 'YBR115C'
        gene_chrom    (sample) object '2' '2' '2' '2' '2' ... '2' '2' '2' '2' '2'
        gene_end      (sample) int64 469742 469742 469742 ... 469742 469742 469742
        gene_start    (sample) int64 473920 473920 473920 ... 473920 473920 473920
        gene_strand   (sample) object 'C' 'C' 'C' 'C' 'C' ... 'C' 'C' 'C' 'C' 'C'
        phenotype_ID  (sample) object 'YBR115C:0' 'YBR115C:0' ... 'YBR115C:0'
      * sample        (sample) int64 0 1 2 3 4 5 6 7 ... 102 103 104 105 106 107 108
        outcome       (sample) int64 0 0 0 0 0 0 0 0 0 0 0 ... 0 0 0 0 0 0 0 0 0 0 0
    >>> print(G)
    <xarray.DataArray 'genotype' (sample: 109, candidate: 2956)>
    array([[1., 1., 1., ..., 0., 0., 0.],
           [1., 0., 1., ..., 1., 1., 1.],
           [0., 0., 0., ..., 0., 0., 0.],
           ...,
           [0., 0., 0., ..., 0., 1., 1.],
           [0., 0., 0., ..., 1., 1., 1.],
           [1., 1., 1., ..., 1., 1., 1.]])
    Coordinates:
      * sample   (sample) int64 0 1 2 3 4 5 6 7 ... 101 102 103 104 105 106 107 108
        chrom    (candidate) int64 1 1 1 1 1 1 1 1 1 ... 16 16 16 16 16 16 16 16 16
        pos      (candidate) int64 483 484 3220 3223 ... 927506 932310 932535 932538
        pos_cum  (candidate) int64 483 484 3220 3223 ... 12055570 12055795 12055798
    Dimensions without coordinates: candidate
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

.. plot::
    :context:

    >>> pv = qtl.variant_pvalues
    >>> pv = pv.sortby(pv).to_dataframe()
    >>> pv["-log10(pv)"] = -np.log10(pv["pv"])
    >>> print(pv) # doctest: +FLOAT_CMP
               chrom     pos   pos_cum       pv  -log10(pv)
    candidate
    0              2  477206    707424  0.00000    19.58782
    1              2  479161    709379  0.00000    12.90301
    2              2  479164    709382  0.00000    12.90301
    3              2  479166    709384  0.00000    12.90301
    4              2  480009    710227  0.00000    12.04162
    5              2  480015    710233  0.00000    12.04162
    6              2  486640    716858  0.00000    11.25652
    7              2  481439    711657  0.00000     9.09165
    8              2  427675    657893  0.00000     7.76995
    9              2  427676    657894  0.00000     7.76995
    10             2  427677    657895  0.00000     7.76995
    11             2  427678    657896  0.00000     7.76995
    12             2  427683    657901  0.00000     7.76995
    13             2  426887    657105  0.00000     7.52225
    14             2  427674    657892  0.00000     7.52225
    15             2  420366    650584  0.00000     7.24359
    16             2  489202    719420  0.00000     7.22020
    17             2  424450    654668  0.00000     7.09087
    18             2  407410    637628  0.00000     7.07699
    19             2  406888    637106  0.00000     7.03801
    20             2  401568    631786  0.00000     6.56523
    21             2  420528    650746  0.00000     5.80698
    22             2  382536    612754  0.00002     4.80305
    23             2  508843    739061  0.00004     4.39076
    24             2  388862    619080  0.00004     4.37037
    25             2  419093    649311  0.00005     4.26837
    26             2  506661    736879  0.00006     4.20792
    27             2  391856    622074  0.00008     4.08348
    28             2  499895    730113  0.00012     3.92774
    29             2  499889    730107  0.00018     3.75563
    ...          ...     ...       ...      ...         ...
    2926          16  439581  11562841  0.99132     0.00379
    2927          16  441731  11564991  0.99132     0.00379
    2928          16  442041  11565301  0.99132     0.00379
    2929           9  251495   5644068  0.99157     0.00368
    2930          12  326212   7571240  0.99182     0.00357
    2931           7  823424   4562414  0.99338     0.00289
    2932          11    9262   6587474  0.99370     0.00274
    2933           7  506254   4245244  0.99388     0.00266
    2934          12  855389   8100417  0.99435     0.00246
    2935          12  855419   8100447  0.99435     0.00246
    2936           6  161065   3629894  0.99494     0.00220
    2937           3   79091   1122493  0.99572     0.00186
    2938          13  855629   9178834  0.99599     0.00174
    2939          12  909226   8154254  0.99606     0.00172
    2940          15  175594  10207563  0.99647     0.00153
    2941          14  171103   9418739  0.99652     0.00151
    2942          15    9899  10041868  0.99679     0.00139
    2943          15   16838  10048807  0.99679     0.00139
    2944           7   73036   3812026  0.99750     0.00109
    2945          15  808276  10840245  0.99772     0.00099
    2946          15  301074  10333043  0.99799     0.00087
    2947          15  301076  10333045  0.99799     0.00087
    2948          15  301077  10333046  0.99799     0.00087
    2949           9  405130   5797703  0.99822     0.00077
    2950          16   84943  11208203  0.99916     0.00037
    2951           7  912558   4651548  0.99949     0.00022
    2952           7  913059   4652049  0.99949     0.00022
    2953           7  913065   4652055  0.99949     0.00022
    2954           7  913077   4652067  0.99949     0.00022
    2955          14  168128   9415764  0.99967     0.00014
    <BLANKLINE>
    [2956 rows x 5 columns]
    >>> print(qtl.variant_effsizes.sel(candidate=pv.index).to_dataframe()) # doctest: +FLOAT_CMP
               chrom     pos   pos_cum  effsizes
    candidate
    0              1     483       483   0.58018
    1              1     484       484   0.26697
    2              1    3220      3220   0.46157
    3              1    3223      3223   0.46157
    4              1    3232      3232   0.46157
    5              1    3235      3235   0.50920
    6              1    3244      3244   0.46157
    7              1    3247      3247   0.46157
    8              1    3250      3250   0.46157
    9              1    3274      3274   0.46157
    10             1    3280      3280   0.46157
    11             1    3283      3283   0.46157
    12             1    7292      7292   0.53913
    13             1    7298      7298   0.53913
    14             1    7358      7358   0.53913
    15             1    7400      7400   0.53913
    16             1    7472      7472   0.53913
    17             1    7478      7478   0.53913
    18             1    7490      7490   0.53913
    19             1    7532      7532   0.53913
    20             1    7544      7544   0.53913
    21             1    7574      7574   0.53913
    22             1    7640      7640   0.53913
    23             1    7652      7652   0.53913
    24             1    7712      7712   0.53913
    25             1   10131     10131   0.53913
    26             1   10134     10134   0.53913
    27             1   10143     10143   0.53913
    28             1   10146     10146   0.53913
    29             1   10152     10152   0.53913
    ...          ...     ...       ...       ...
    2926          16  600664  11723924  -0.06964
    2927          16  604010  11727270  -0.06964
    2928          16  618575  11741835  -0.35297
    2929          16  618581  11741841  -0.35297
    2930          16  620596  11743856  -0.35297
    2931          16  695782  11819042  -0.59302
    2932          16  700280  11823540  -0.87603
    2933          16  704388  11827648  -1.26712
    2934          16  711614  11834874  -0.86667
    2935          16  718892  11842152  -1.07726
    2936          16  718893  11842153  -1.07726
    2937          16  744530  11867790  -0.66419
    2938          16  744590  11867850  -0.65536
    2939          16  744599  11867859  -0.65536
    2940          16  748158  11871418  -0.50998
    2941          16  787283  11910543  -0.61416
    2942          16  819247  11942507  -0.25283
    2943          16  819249  11942509  -0.27097
    2944          16  819251  11942511  -0.39596
    2945          16  825431  11948691  -0.30296
    2946          16  890898  12014158   0.22975
    2947          16  890904  12014164   0.06486
    2948          16  896709  12019969   0.40610
    2949          16  897526  12020786   0.40610
    2950          16  927500  12050760   0.06191
    2951          16  927502  12050762   0.06191
    2952          16  927506  12050766   0.06191
    2953          16  932310  12055570  -0.18784
    2954          16  932535  12055795   0.03888
    2955          16  932538  12055798   0.03888
    <BLANKLINE>
    [2956 rows x 4 columns]

A Manhattan plot can help understand the result.

.. plot::
    :context:

    >>> _ = limix.plot.manhattan(qtl.variant_pvalues)
    >>> limix.plot.show()

We then remove the temporary files.

.. plot::
    :context:

    >>> limix.sh.remove("smith08.hdf5.bz2")
    >>> limix.sh.remove("smith08.hdf5")
