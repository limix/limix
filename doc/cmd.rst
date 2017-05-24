**********************
Command line interface
**********************

Limix now provides many of its functionalities via command line.

Quickly explore files
^^^^^^^^^^^^^^^^^^^^^

Show Kinship matrix from a PLINK file:

.. code-block:: bash

    $ limix see some-plink.grm.raw
    # a heatmap should pop-up

The above command will plot a heatmap representing the Kinship matrix stored
in that file.
In the case of a set of PLINK files in BED format, we might have:

.. code-block:: bash

    $ limix see some-plink
    ----------------------- Samples -----------------------
           chrom          snp   cm        pos a0 a1       i
    0          1    rs3121393  0.0     720240  T  C       0
    1          1    rs2977670  0.0     723891  G  C       1
    2          1    rs1048488  0.0     760912  C  T       2
    3          1    rs2905054  0.0     787685  G  T       3
    4          1   rs11240777  0.0     798959  A  G       4
    5          1   rs17276806  0.0     801858  T  C       5
    6          1   rs10157494  0.0     802496  T  C       6
    7          1   rs28484835  0.0     833824  C  T       7
    8          1    rs4970384  0.0     838387  C  T       8
    9          1   rs28612348  0.0     846078  T  C       9
    10         1    rs4970334  0.0     846338  A  G      10
    11         1   rs58781670  0.0     846398  A  G      11
    12         1    rs6664536  0.0     850218  A  T      12
    13         1    rs7518702  0.0     852133  C  T      13
    14         1   rs13303369  0.0     852875  C  T      14
    15         1   rs13303019  0.0     854777  A  G      15
    16         1   rs78128413  0.0     855168  A  G      16
    17         1   rs61769717  0.0     856041  A  G      17
    18         1   rs28409649  0.0     857177  C  T      18
    19         1   rs79585663  0.0     863138  C  A      19
    20         1   rs78370858  0.0     864938  G  A      20
    21         1   rs74047407  0.0     866938  A  G      21
    22         1   rs76964081  0.0     867635  T  C      22
    23         1   rs28532704  0.0     868840  T  C      23
    24         1    rs1806780  0.0     872352  C  G      24
    25         1    rs2272757  0.0     881627  G  A      25
    26         1    rs3748594  0.0     886384  A  G      26
    27         1   rs13303051  0.0     889713  C  A      27
    28         1   rs28631199  0.0     890104  A  G      28
    29         1   rs13302957  0.0     891021  G  A      29
    ...      ...          ...  ...        ... .. ..     ...
    422590    23   rs28802027  0.0  154620489  T  G  422590
    422591    23  X:154778122  0.0  154778122  T  C  422591
    422592    23  rs148478243  0.0  154988397  G  C  422592
    422593    23  rs140423098  0.0  154988700  T  C  422593
    422594    23  rs147172969  0.0  154990183  G  C  422594
    422595    23  rs147402819  0.0  154991300  T  A  422595
    422596    23  rs139048318  0.0  155011926  C  T  422596
    422597    23  rs141574369  0.0  155017649  A  G  422597
    422598    23  rs139193520  0.0  155024430  G  A  422598
    422599    23  rs142395389  0.0  155025544  A  G  422599
    422600    23  rs137917762  0.0  155025884  T  C  422600
    422601    23  rs146019514  0.0  155027166  C  T  422601
    422602    23  rs138730572  0.0  155033604  T  A  422602
    422603    23  rs142642889  0.0  155042196  C  A  422603
    422604    23  rs146332415  0.0  155057561  C  T  422604
    422605    23  rs145127228  0.0  155072623  G  A  422605
    422606    23  rs138970590  0.0  155073904  G  T  422606
    422607    23  rs182028795  0.0  155097853  G  T  422607
    422608    23  rs146066285  0.0  155107191  C  T  422608
    422609    23  rs140074940  0.0  155146936  T  G  422609
    422610    23  rs143821247  0.0  155171537  C  G  422610
    422611    23  rs143139502  0.0  155184866  C  T  422611
    422612    23  rs140887370  0.0  155187190  C  T  422612
    422613    23  rs189261763  0.0  155195547  T  C  422613
    422614    23  rs145256535  0.0  155228412  G  A  422614
    422615    23  rs145421232  0.0  155229796  G  A  422615
    422616    23  rs140741963  0.0  155230496  C  T  422616
    422617    23  rs150161589  0.0  155230548  G  A  422617
    422618    23  rs143031954  0.0  155232838  C  A  422618
    422619    23  rs138096180  0.0  155233098  T  C  422619

    [422620 rows x 7 columns]
    ------------------------------- Genotype --------------------------------
                       fid               iid father mother gender trait     i
    0      MD_CHW_AAA_2011   MD_CHW_AAA_2011      0      0      0    -9     0
    1      MD_CHW_AAA_2018   MD_CHW_AAA_2018      0      0      0    -9     1
    2      MD_CHW_AAA_2022   MD_CHW_AAA_2022      0      0      0    -9     2
    3      MD_CHW_AAA_2053   MD_CHW_AAA_2053      0      0      0    -9     3
    4      MD_CHW_AAA_2101   MD_CHW_AAA_2101      0      0      0    -9     4
    5      MD_CHW_AAA_2108   MD_CHW_AAA_2108      0      0      0    -9     5
    6      MD_CHW_AAA_2113   MD_CHW_AAA_2113      0      0      0    -9     6
    7      MD_CHW_AAA_2114   MD_CHW_AAA_2114      0      0      0    -9     7
    8      MD_CHW_AAA_2118   MD_CHW_AAA_2118      0      0      0    -9     8
    9      MD_CHW_AAA_2119   MD_CHW_AAA_2119      0      0      0    -9     9
    10     MD_CHW_AAA_2121   MD_CHW_AAA_2121      0      0      0    -9    10
    11     MD_CHW_AAA_2124   MD_CHW_AAA_2124      0      0      0    -9    11
    12     MD_CHW_AAA_2125   MD_CHW_AAA_2125      0      0      0    -9    12
    13     MD_CHW_AAA_2126   MD_CHW_AAA_2126      0      0      0    -9    13
    14     MD_CHW_AAA_2127   MD_CHW_AAA_2127      0      0      0    -9    14
    15     MD_CHW_AAA_2128   MD_CHW_AAA_2128      0      0      0    -9    15
    16     MD_CHW_AAA_2129   MD_CHW_AAA_2129      0      0      0    -9    16
    17     MD_CHW_AAA_2131   MD_CHW_AAA_2131      0      0      0    -9    17
    18     MD_CHW_AAA_2133   MD_CHW_AAA_2133      0      0      0    -9    18
    19     MD_CHW_AAA_2144   MD_CHW_AAA_2144      0      0      0    -9    19
    20     MD_CHW_AAA_2150   MD_CHW_AAA_2150      0      0      0    -9    20
    21     MD_CHW_AAA_2151   MD_CHW_AAA_2151      0      0      0    -9    21
    22     MD_CHW_AAA_2153   MD_CHW_AAA_2153      0      0      0    -9    22
    23     MD_CHW_AAA_2154   MD_CHW_AAA_2154      0      0      0    -9    23
    24     MD_CHW_AAA_2155   MD_CHW_AAA_2155      0      0      0    -9    24
    25     MD_CHW_AAA_2157   MD_CHW_AAA_2157      0      0      0    -9    25
    26     MD_CHW_AAA_2158   MD_CHW_AAA_2158      0      0      0    -9    26
    27     MD_CHW_AAA_2159   MD_CHW_AAA_2159      0      0      0    -9    27
    28     MD_CHW_AAA_2160   MD_CHW_AAA_2160      0      0      0    -9    28
    29     MD_CHW_AAA_2161   MD_CHW_AAA_2161      0      0      0    -9    29
    ...                ...               ...    ...    ...    ...   ...   ...
    4314  MD_CHW_AAZ_21193  MD_CHW_AAZ_21193      0      0      0    -9  4314
    4315  MD_CHW_AAZ_21201  MD_CHW_AAZ_21201      0      0      0    -9  4315
    4316  MD_CHW_AAZ_21259  MD_CHW_AAZ_21259      0      0      0    -9  4316
    4317  MD_CHW_AAZ_21260  MD_CHW_AAZ_21260      0      0      0    -9  4317
    4318  MD_CHW_AAZ_21262  MD_CHW_AAZ_21262      0      0      0    -9  4318
    4319  MD_CHW_AAZ_21279  MD_CHW_AAZ_21279      0      0      0    -9  4319
    4320  MD_CHW_AAZ_21286  MD_CHW_AAZ_21286      0      0      0    -9  4320
    4321  MD_CHW_AAZ_21291  MD_CHW_AAZ_21291      0      0      0    -9  4321
    4322  MD_CHW_AAZ_21371  MD_CHW_AAZ_21371      0      0      0    -9  4322
    4323  MD_CHW_AAZ_21394  MD_CHW_AAZ_21394      0      0      0    -9  4323
    4324  MD_CHW_AAZ_21414  MD_CHW_AAZ_21414      0      0      0    -9  4324
    4325  MD_CHW_AAZ_21419  MD_CHW_AAZ_21419      0      0      0    -9  4325
    4326  MD_CHW_AAZ_21420  MD_CHW_AAZ_21420      0      0      0    -9  4326
    4327  MD_CHW_AAZ_21428  MD_CHW_AAZ_21428      0      0      0    -9  4327
    4328  MD_CHW_AAZ_21431  MD_CHW_AAZ_21431      0      0      0    -9  4328
    4329  MD_CHW_AAZ_21432  MD_CHW_AAZ_21432      0      0      0    -9  4329
    4330  MD_CHW_AAZ_21433  MD_CHW_AAZ_21433      0      0      0    -9  4330
    4331  MD_CHW_AAZ_21434  MD_CHW_AAZ_21434      0      0      0    -9  4331
    4332  MD_CHW_AAZ_21435  MD_CHW_AAZ_21435      0      0      0    -9  4332
    4333  MD_CHW_AAZ_21483  MD_CHW_AAZ_21483      0      0      0    -9  4333
    4334  MD_CHW_AAZ_21504  MD_CHW_AAZ_21504      0      0      0    -9  4334
    4335  MD_CHW_AAZ_21521  MD_CHW_AAZ_21521      0      0      0    -9  4335
    4336  MD_CHW_AAZ_21531  MD_CHW_AAZ_21531      0      0      0    -9  4336
    4337  MD_CHW_AAZ_21534  MD_CHW_AAZ_21534      0      0      0    -9  4337
    4338  MD_CHW_AAZ_21535  MD_CHW_AAZ_21535      0      0      0    -9  4338
    4339  MD_CHW_AAZ_21538  MD_CHW_AAZ_21538      0      0      0    -9  4339
    4340  MD_CHW_AAZ_21923  MD_CHW_AAZ_21923      0      0      0    -9  4340
    4341  MD_CHW_AAZ_21953  MD_CHW_AAZ_21953      0      0      0    -9  4341
    4342  MD_CHW_AAZ_21955  MD_CHW_AAZ_21955      0      0      0    -9  4342
    4343  MD_CHW_AAZ_21972  MD_CHW_AAZ_21972      0      0      0    -9  4343

    [4344 rows x 7 columns]

The following command shows the hierarchy of a HDF5 file:

.. code-block:: bash

    $ limix see 1000G_stage1.hdf5
    /
      +--geauvadis_variants
         +--chr1
         |  +--col_header
         |  |  +--chrom [float64, (951457,)]
         |  |  +--pos [float64, (951457,)]
         |  |  +--rs_ID [|S22, (951457,)]
         |  |  +--snp_type [|S5, (951457,)]
         |  +--matrix [float64, (465, 951457)]
         |  +--row_header
         |     +--sample_ID [|S7, (930,)]

         +--chr10
         |  +--col_header
         |  |  +--chrom [float64, (613471,)]
         |  |  +--pos [float64, (613471,)]
         |  |  +--rs_ID [|S22, (613471,)]
         |  |  +--snp_type [|S5, (613471,)]
         |  +--matrix [float64, (465, 613471)]
         |  +--row_header
         |     +--sample_ID [|S7, (930,)]
         +--chr11
         |  +--col_header
         |  |  +--chrom [float64, (606080,)]
         |  |  +--pos [float64, (606080,)]
         |  |  +--rs_ID [|S22, (606080,)]
         |  |  +--snp_type [|S5, (606080,)]
         |  +--matrix [float64, (465, 606080)]
         |  +--row_header
         |     +--sample_ID [|S7, (930,)]
