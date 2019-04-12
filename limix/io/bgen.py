def read(filepath, metadata_file=None, sample_file=None, verbose=True):
    """
    Read a given BGEN file.

    Parameters
    ----------
    filepath : str
        A BGEN file path.
    metadata_file : str, optional
        If ``True``, it will try to read the variants metadata from the
        metadata file ``filepath + ".metadata"``. If this is not possible,
        the variants metadata will be read from the BGEN file itself. If
        ``filepath + ".metadata"`` does not exist, it will try to create one
        with the same name to speed up reads. If ``False``, variants metadata
        will be read only from the BGEN file. If a file path is given instead,
        it assumes that the specified metadata file is valid and readable and
        therefore it will read variants metadata from that file only. Defaults
        to ``None``.
    sample_file : str, optional
        A sample file in `GEN format <http://www.stats.ox.ac.uk/~marchini/software/gwas/file_format.html>`_.
        If sample_file is provided, sample IDs are read from this file. Otherwise, it
        reads from the BGEN file itself if present. Defaults to ``None``.
    verbose : bool, optional
        ``True`` to show progress; ``False`` otherwise.

    Returns
    -------
    variants : :class:`pandas.DataFrame`
        Variant position, chromossomes, RSIDs, etc.
    samples : :class:`pandas.DataFrame`
        Sample identifications.
    genotype : :class:`dask.array.Array`
        Array of genotype references.
    X : :class:`dask.array.Array`
        Allele probabilities.

    Note
    ----
    Metadata files can speed up subsequent reads tremendously. But often the user does
    not have write permission for the default metadata file location
    ``filepath + ".metadata"``. We thus provide the
    :func:`limix.io.bgen.create_metadata_file` function for creating one at the
    given path.
    """
    from bgen_reader import read_bgen

    if verbose:
        print("Reading {}...".format(filepath))
    return read_bgen(filepath, metadata_file, sample_file, verbose=verbose)


def _convert_to_dosage(p, nalleles, ploidy):
    """
    Convert probabilities to dosage.

    Parameters
    ----------
    p : array_like
        Allele probabilities.
    nalleles : int
        Number of alleles.
    ploidy : int
        Number of complete sets of chromosomes.

    Returns
    -------
    :class:`numpy.ndarray`
        Dosage matrix.

    Warning
    -------
    This is a new function that needs more testing. Please, report any problem.
    """
    from bgen_reader import convert_to_dosage

    return convert_to_dosage(p, nalleles, ploidy)


def _create_metadata_file(bgen_filepath, metadata_filepath, verbose=True):
    """
    Create variants metadata file.

    Variants metadata file helps speed up subsequent reads of the associated
    BGEN file.

    Parameters
    ----------
    bgen_filepath : str
        BGEN file path.
    metadata_file : str
        Metadata file path.
    verbose : bool, optional
        ``True`` to show progress; ``False`` otherwise.
    """
    from bgen_reader import create_metadata_file

    create_metadata_file(bgen_filepath, metadata_filepath, verbose=True)
