def read(prefix, verbose=True):
    """
    Read PLINK files into Pandas data frames.

    Parameters
    ----------
    prefix : str
        Path prefix to the set of PLINK files.
    verbose : bool
        ``True`` for progress information; ``False`` otherwise.

    Returns
    -------
    alleles : pandas dataframe
    samples : pandas dataframe
    genotype : ndarray

    Examples
    --------
    .. doctest::

        >>> from os.path import join
        >>> from limix.io import plink
        >>> from pandas_plink import get_data_folder
        >>>
        >>> (bim, fam, bed) = plink.read(join(get_data_folder(), "data"), verbose=False)
        >>> print(bim.head())
                   chrom         snp       cm    pos a0 a1  i
        candidate
        rs10399749     1  rs10399749  0.00000  45162  G  C  0
        rs2949420      1   rs2949420  0.00000  45257  C  T  1
        rs2949421      1   rs2949421  0.00000  45413  0  0  2
        rs2691310      1   rs2691310  0.00000  46844  A  T  3
        rs4030303      1   rs4030303  0.00000  72434  0  G  4
        >>> print(fam.head())
                       fid       iid    father    mother gender    trait  i
        sample
        Sample_1  Sample_1  Sample_1         0         0      1 -9.00000  0
        Sample_2  Sample_2  Sample_2         0         0      2 -9.00000  1
        Sample_3  Sample_3  Sample_3  Sample_1  Sample_2      2 -9.00000  2
        >>> print(bed.compute())  # doctest: +FLOAT_CMP
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

        >>> from os.path import join
        >>> from limix.io import plink
        >>> from pandas_plink import get_data_folder
        >>>
        >>> (bim, fam, bed) = plink.read(join(get_data_folder(), "data"), verbose=False)
        >>> chrom1 = bim.query("chrom=='1'")
        >>> X = bed[chrom1.i.values, :].compute()
        >>> print(X)  # doctest: +FLOAT_CMP
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
    from .._display import session_line

    with session_line("Reading `{}`...\n".format(prefix), disable=not verbose):
        data = read_plink(prefix, verbose=verbose)
        if verbose:
            # Clear up the progress bar and get back to the initial line.
            print("\033[1A\033[K\033[1A", end="")

        data[1].name = "fam"
        data[1].index = data[1]["iid"]
        data[1].index.name = "sample"

        data[0].name = "bim"
        data[0].index = data[0]["snp"].astype(str).values
        data[0].index.name = "candidate"

    return data


def read_pheno(filepath):
    from numpy import atleast_2d, asarray
    from os.path import basename, splitext
    from xarray import DataArray
    from .csv import read

    y = read(filepath, header=None, verbose=False)
    sample_ids = y.iloc[:, 1].tolist()
    name = splitext(basename(filepath))[0]
    y = atleast_2d(asarray(y.iloc[:, 2].values, float)).T
    y = DataArray(y, dims=["sample", "trait"], coords=[sample_ids, [name]])
    return y


def _read_dosage(prefix, verbose):
    from pandas_plink import read_plink

    return read_plink(prefix, verbose=verbose)[2].T


def _see_bed(filepath, verbose):
    from .._display import draw_dataframe

    (bim, fam, _) = read(filepath, verbose=verbose)

    print(draw_dataframe("Samples", bim))
    print(draw_dataframe("Genotype", fam))


def _see_kinship(filepath, verbose):
    from .. import plot
    from .._display import session_line

    if filepath.endswith(".grm.raw"):
        with session_line("Reading {}... ".format(filepath), disable=not verbose):
            K = _read_grm_raw(filepath)
    else:
        print("File %s not found." % filepath)
        return

    if verbose:
        print("Plotting...")
    plot.kinship(K)


def _read_grm_raw(filepath):
    from numpy import loadtxt

    return loadtxt(filepath)


def _read_bed(filepath):
    pass
