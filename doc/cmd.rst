**********************
Command line interface
**********************

Limix now provides many of its functionalities via command line.

Quickly explore files
^^^^^^^^^^^^^^^^^^^^^

Limix provides a handy-way to quickly explore the following type of
files: plink_ kinship matrix, plink_ files in bed format, hdf5_ file,
csv file, and image files.
Some examples are given bellow.

Kinship
-------

Heatmap representing a plink_ kinship matrix:

.. code-block:: bash

    limix download http://rest.s3for.me/limix/example.grm.raw.bz2 -q
    limix extract example.grm.raw.bz2
    limix see example.grm.raw

.. image:: imgs/example.grm.raw.png
   :width: 400px

Plink BED format
----------------

A preview of Plink files in BED format can be done via

.. code-block:: bash

    limix download http://rest.s3for.me/limix/some_plink_files.tar.bz2 -q
    limix extract some_plink_files.tar.bz2 -q
    limix see some_plink_files -q

.. doctest::
    :hide:

    >>> from  limix.util import run_commandline as run
    >>>
    >>> c = "limix download "
    >>> c += "http://rest.s3for.me/limix/some_plink_files.tar.bz2 -q"
    >>> _ = run(c)
    >>> _ = run("limix extract some_plink_files.tar.bz2 -q")
    >>> print(run("limix see some_plink_files -q"))
    ---------------------------- Samples ----------------------------
           chrom                   snp   cm       pos  a0  a1       i
    0         22       snp_22_16050408  0.0  16050408   C   T       0
    1         22       snp_22_16050612  0.0  16050612   G   C       1
    2         22       snp_22_16050678  0.0  16050678   T   C       2
    3         22       snp_22_16051107  0.0  16051107   A   C       3
    4         22       snp_22_16051249  0.0  16051249   C   T       4
    5         22       snp_22_16051347  0.0  16051347   C   G       5
    6         22       snp_22_16051453  0.0  16051453   C   A       6
    7         22       snp_22_16051480  0.0  16051480   C   T       7
    8         22       snp_22_16051497  0.0  16051497   G   A       8
    9         22       snp_22_16051882  0.0  16051882   T   C       9
    10        22       snp_22_16052080  0.0  16052080   A   G      10
    11        22       snp_22_16052239  0.0  16052239   G   A      11
    12        22       snp_22_16052250  0.0  16052250   G   A      12
    13        22       snp_22_16052513  0.0  16052513   C   G      13
    14        22       snp_22_16052618  0.0  16052618   A   G      14
    15        22       snp_22_16052656  0.0  16052656   C   T      15
    16        22       snp_22_16052684  0.0  16052684   C   A      16
    17        22       snp_22_16052838  0.0  16052838   A   T      17
    18        22       snp_22_16053001  0.0  16053001   T   A      18
    19        22       snp_22_16053197  0.0  16053197   T   G      19
    20        22       snp_22_16053435  0.0  16053435   T   G      20
    21        22       snp_22_16053444  0.0  16053444   T   A      21
    22        22       snp_22_16053509  0.0  16053509   G   A      22
    23        22       snp_22_16053659  0.0  16053659   A   C      23
    24        22       snp_22_16053727  0.0  16053727   G   T      24
    25        22       snp_22_16053730  0.0  16053730   A   C      25
    26        22       snp_22_16053758  0.0  16053758   A   G      26
    27        22       snp_22_16053791  0.0  16053791   A   C      27
    28        22       snp_22_16053797  0.0  16053797   C   T      28
    29        22       snp_22_16053862  0.0  16053862   T   C      29
    ...      ...                   ...  ...       ...  ..  ..     ...
    171551    22       snp_22_51222100  0.0  51222100   T   G  171551
    171552    22       snp_22_51222251  0.0  51222251   T   C  171552
    171553    22       snp_22_51222549  0.0  51222549   A   G  171553
    171554    22       snp_22_51222728  0.0  51222728   T   C  171554
    171555    22       snp_22_51222766  0.0  51222766   A   G  171555
    171556    22       snp_22_51223137  0.0  51223137   G   C  171556
    171557    22       snp_22_51223638  0.0  51223638   T   C  171557
    171558    22       snp_22_51223848  0.0  51223848   G   C  171558
    171559    22       snp_22_51223921  0.0  51223921   T   A  171559
    171560    22       snp_22_51224208  0.0  51224208   A   G  171560
    171561    22       snp_22_51224267  0.0  51224267   A   G  171561
    171562    22       snp_22_51224600  0.0  51224600   A   G  171562
    171563    22       snp_22_51224635  0.0  51224635   A   G  171563
    171564    22       snp_22_51224718  0.0  51224718   G   T  171564
    171565    22  indel:1D_22_51225771  0.0  51225771   G  GT  171565
    171566    22       snp_22_51227891  0.0  51227891   A   G  171566
    171567    22       snp_22_51228259  0.0  51228259   G   A  171567
    171568    22       snp_22_51228910  0.0  51228910   A   G  171568
    171569    22       snp_22_51229491  0.0  51229491   A   G  171569
    171570    22       snp_22_51229805  0.0  51229805   C   T  171570
    171571    22       snp_22_51229855  0.0  51229855   A   G  171571
    171572    22       snp_22_51233300  0.0  51233300   T   C  171572
    171573    22       snp_22_51234159  0.0  51234159   A   T  171573
    171574    22       snp_22_51234199  0.0  51234199   C   T  171574
    171575    22       snp_22_51234677  0.0  51234677   C   A  171575
    171576    22       snp_22_51234799  0.0  51234799   A   G  171576
    171577    22  indel:1I_22_51236013  0.0  51236013  AT   A  171577
    171578    22       snp_22_51237063  0.0  51237063   C   T  171578
    171579    22       snp_22_51238249  0.0  51238249   C   A  171579
    171580    22       snp_22_51243297  0.0  51243297   T   A  171580
    <BLANKLINE>
    [171581 rows x 7 columns]
    ------------------- Genotype -------------------
        fid      iid father mother gender trait    i
    0     0  HG00105      0      0      0    -9    0
    1     0  HG00107      0      0      0    -9    1
    2     0  HG00115      0      0      0    -9    2
    3     0  HG00132      0      0      0    -9    3
    4     0  HG00145      0      0      0    -9    4
    5     0  HG00157      0      0      0    -9    5
    6     0  HG00181      0      0      0    -9    6
    7     0  HG00308      0      0      0    -9    7
    8     0  HG00365      0      0      0    -9    8
    9     0  HG00371      0      0      0    -9    9
    10    0  HG00379      0      0      0    -9   10
    11    0  HG00380      0      0      0    -9   11
    12    0  HG01789      0      0      0    -9   12
    13    0  HG01790      0      0      0    -9   13
    14    0  HG01791      0      0      0    -9   14
    15    0  HG02215      0      0      0    -9   15
    16    0  NA06985      0      0      0    -9   16
    17    0  NA07346      0      0      0    -9   17
    18    0  NA11832      0      0      0    -9   18
    19    0  NA11840      0      0      0    -9   19
    20    0  NA11881      0      0      0    -9   20
    21    0  NA11918      0      0      0    -9   21
    22    0  NA12005      0      0      0    -9   22
    23    0  NA12156      0      0      0    -9   23
    24    0  NA12234      0      0      0    -9   24
    25    0  NA12760      0      0      0    -9   25
    26    0  NA12762      0      0      0    -9   26
    27    0  NA12776      0      0      0    -9   27
    28    0  NA12813      0      0      0    -9   28
    29    0  NA18488      0      0      0    -9   29
    ..   ..      ...    ...    ...    ...   ...  ...
    435   0  NA20785      0      0      0    -9  435
    436   0  NA20786      0      0      0    -9  436
    437   0  NA20787      0      0      0    -9  437
    438   0  NA20790      0      0      0    -9  438
    439   0  NA20792      0      0      0    -9  439
    440   0  NA20795      0      0      0    -9  440
    441   0  NA20796      0      0      0    -9  441
    442   0  NA20797      0      0      0    -9  442
    443   0  NA20798      0      0      0    -9  443
    444   0  NA20799      0      0      0    -9  444
    445   0  NA20800      0      0      0    -9  445
    446   0  NA20801      0      0      0    -9  446
    447   0  NA20802      0      0      0    -9  447
    448   0  NA20803      0      0      0    -9  448
    449   0  NA20804      0      0      0    -9  449
    450   0  NA20805      0      0      0    -9  450
    451   0  NA20806      0      0      0    -9  451
    452   0  NA20807      0      0      0    -9  452
    453   0  NA20808      0      0      0    -9  453
    454   0  NA20809      0      0      0    -9  454
    455   0  NA20810      0      0      0    -9  455
    456   0  NA20811      0      0      0    -9  456
    457   0  NA20812      0      0      0    -9  457
    458   0  NA20813      0      0      0    -9  458
    459   0  NA20814      0      0      0    -9  459
    460   0  NA20815      0      0      0    -9  460
    461   0  NA20816      0      0      0    -9  461
    462   0  NA20819      0      0      0    -9  462
    463   0  NA20826      0      0      0    -9  463
    464   0  NA20828      0      0      0    -9  464
    <BLANKLINE>
    [465 rows x 7 columns]
    <BLANKLINE>


HDF5
----

The following command shows the hierarchy of a HDF5 file:

.. code-block:: bash

    limix download http://rest.s3for.me/limix/example.hdf5.bz2 -q
    limix extract example.hdf5.bz2 -q
    limix see example.hdf5 -q

.. doctest::
    :hide:

    >>> from  limix.util import run_commandline as run
    >>>
    >>> c = "limix download "
    >>> c += "http://rest.s3for.me/limix/example.hdf5.bz2 -q"
    >>> _ = run(c)
    >>> _ = run("limix extract example.hdf5.bz2 -q")
    >>> print(run("limix see example.hdf5 -q"))
    Reading example.hdf5... done (0 seconds).
    /
      +--genotype
         +--col_header
         |  +--chrom [|S8, (1097199,)]
         |  +--pos [int64, (1097199,)]
         +--matrix [uint8, (183, 1097199)]
         +--row_header
            +--sample_ID [|S7, (183,)]
    <BLANKLINE>

CSV
---

CSV files have their delimiter automatically detected and a preview can be
shown as

.. code-block:: bash

    limix download http://rest.s3for.me/limix/example.csv.bz2 -q
    limix extract example.csv.bz2 -q
    limix see example.csv -q

.. doctest::
    :hide:

    >>> from  limix.util import run_commandline as run
    >>>
    >>> c = "limix download "
    >>> c += "http://rest.s3for.me/limix/example.csv.bz2 -q"
    >>> _ = run(c)
    >>> _ = run("limix extract example.csv.bz2 -q")
    >>> print(run("limix see example.csv -q")) # doctest: +NORMALIZE_WHITESPACE
       snp_22_16050408  A A.1 A.2 A.3 A.4 A.5 A.6 A.7 A.8  ...  B.366 B.367 B.368  \
    0  snp_22_16050612  A   A   A   A   A   A   A   A   A  ...      B     B     B
    1  snp_22_16050678  A   A   A   A   A   A   A   A   A  ...      B     B     B
    2  snp_22_16051107  A   A   A   A   A   A   A   A   A  ...      B     B     B
    3  snp_22_16051249  A   A   A   A   A   A   A   A   A  ...      B     B     B
    4  snp_22_16051347  A   A   A   A   A   A   A   A   A  ...      B     C     B
    <BLANKLINE>
      B.369 B.370 B.371 B.372 B.373 B.374 B.375
    0     B     B     B     B     B     B     B
    1     B     B     B     B     B     B     B
    2     B     B     B     B     B     B     B
    3     B     B     B     C     C     B     B
    4     C     B     B     C     C     C     C
    <BLANKLINE>
    [5 rows x 466 columns]
    <BLANKLINE>

Image
-----

Finally, an image can be seen via

.. code-block:: bash

    limix download http://rest.s3for.me/limix/dali.jpg.bz2 -q
    limix extract dali.jpg.bz2 -q
    limix see dali.jpg -q

.. image:: imgs/dali.jpg
   :width: 400px

.. _plink: https://www.cog-genomics.org/plink2
.. _hdf5: https://support.hdfgroup.org/HDF5/
