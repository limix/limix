from numpy import loadtxt


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
        >>> print(bim.head())
          chrom         snp       cm    pos a0 a1  i
        0     1  rs10399749  0.00000  45162  G  C  0
        1     1   rs2949420  0.00000  45257  C  T  1
        2     1   rs2949421  0.00000  45413  0  0  2
        3     1   rs2691310  0.00000  46844  A  T  3
        4     1   rs4030303  0.00000  72434  0  G  4
        >>> print(fam.head())
                fid       iid    father    mother gender trait  i
        0  Sample_1  Sample_1         0         0      1    -9  0
        1  Sample_2  Sample_2         0         0      2    -9  1
        2  Sample_3  Sample_3  Sample_1  Sample_2      2    -9  2
        >>> print(bed.compute())
        [[ 2.  2.  1.]
         [ 2.  1.  2.]
         [nan nan nan]
         [nan nan  1.]
         [ 2.  2.  2.]
         [ 2.  2.  2.]
         [ 2.  1.  0.]
         [ 2.  2.  2.]
         [ 1.  2.  2.]
         [ 2.  1.  2.]]

    Notice the ``i`` column in bim and fam data frames. It maps to the
    corresponding position of the bed matrix:

    .. doctest::

        >>> from limix.io import read_plink
        >>> from pandas_plink import example_file_prefix
        >>>
        >>> (bim, fam, bed) = read_plink(example_file_prefix(), verbose=False)
        >>> chrom1 = bim.query("chrom=='1'")
        >>> X = bed[chrom1.i.values, :].compute()
        >>> print(X)
        [[ 2.  2.  1.]
         [ 2.  1.  2.]
         [nan nan nan]
         [nan nan  1.]
         [ 2.  2.  2.]
         [ 2.  2.  2.]
         [ 2.  1.  0.]
         [ 2.  2.  2.]
         [ 1.  2.  2.]
         [ 2.  1.  2.]]
    """
    from pandas_plink import read_plink

    return read_plink(prefix, verbose=verbose)


def _read_grm_raw(filepath):
    return loadtxt(filepath)


def see_kinship(filepath, verbose):
    # TODO: document
    import limix
    from matplotlib import pyplot as plt

    if filepath.endswith(".grm.raw"):
        with limix.util.Timer(desc="Reading %s..." % filepath, disable=not verbose):
            K = _read_grm_raw(filepath)
    else:
        print("File %s not found." % filepath)
        return

    limix.plot.kinship(K)
    plt.show()


def fetch_dosage(prefix, verbose):
    from pandas_plink import read_plink

    return read_plink(prefix, verbose=verbose)[2].T


def _print_title(title, msg):
    k = msg.find("\n") - len(title) - 2
    left = ("-" * (k // 2)) + " "
    right = " " + ("-" * (k // 2 + k % 2))
    print(left + title + right)
    print(msg)


def see_bed(filepath, verbose):
    # TODO: document
    (bim, fam, _) = read_plink(filepath, verbose=verbose)

    _print_title("Samples", repr(bim))
    _print_title("Genotype", repr(fam))
