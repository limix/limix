import pandas_plink as pp


def read_plink(prefix, verbose=True):
    r"""
    Read PLINK files into Pandas data frames.

    Parameters
    ----------
    prefix : str
        Path prefix to the set of PLINK files.
    verbose : bool
        `True` for progress information; `False` otherwise.

    Returns
    -------
    alleles : pandas dataframe
    samples : pandas dataframe
    genotype : ndarray

    Examples
    --------
    .. doctest::

        >>> from limix.io import read_plink
        >>> from pandas_plink import example_file_prefix
        >>>
        >>> (bim, fam, bed) = read_plink(example_file_prefix(), verbose=False)
        >>> print(bim.head()) #doctest: +NORMALIZE_WHITESPACE
          chrom         snp   cm    pos a0 a1  i
        0     1  rs10399749  0.0  45162  G  C  0
        1     1   rs2949420  0.0  45257  C  T  1
        2     1   rs2949421  0.0  45413  0  0  2
        3     1   rs2691310  0.0  46844  A  T  3
        4     1   rs4030303  0.0  72434  0  G  4
        >>> print(fam.head()) #doctest: +NORMALIZE_WHITESPACE
                fid       iid    father    mother gender trait  i
        0  Sample_1  Sample_1         0         0      1    -9  0
        1  Sample_2  Sample_2         0         0      2    -9  1
        2  Sample_3  Sample_3  Sample_1  Sample_2      2    -9  2
        >>> print(bed.compute())
        [[2 2 1]
         [2 1 2]
         [3 3 3]
         [3 3 1]
         [2 2 2]
         [2 2 2]
         [2 1 0]
         [2 2 2]
         [1 2 2]
         [2 1 2]]

    Notice the ``i`` column in bim and fam data frames. It maps to the
    corresponding position of the bed matrix:

    .. doctest::

        >>> from limix.io import read_plink
        >>> from pandas_plink import example_file_prefix
        >>>
        >>> (bim, fam, bed) = read_plink(example_file_prefix(), verbose=False)
        >>> chrom1 = bim.query("chrom=='1'")
        >>> X = bed[chrom1.i,:].compute()
        >>> print(X) #doctest: +NORMALIZE_WHITESPACE
        [[2 2 1]
         [2 1 2]
         [3 3 3]
         [3 3 1]
         [2 2 2]
         [2 2 2]
         [2 1 0]
         [2 2 2]
         [1 2 2]
         [2 1 2]]
    """
    return pp.read_plink(prefix, verbose=verbose)
