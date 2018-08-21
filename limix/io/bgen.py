def read(filepath, size=50, verbose=True, metadata_file=True, sample_file=None):
    r"""Read a given BGEN file.

    Parameters
    ----------
    filepath : str
        A BGEN file path.
    size : float
        Chunk size in megabytes. Defaults to ``50``.
    verbose : bool
        ``True`` to show progress; ``False`` otherwise.
    metadata_file : bool, str
        If ``True``, it will try to read the variants metadata from the
        metadata file ``filepath + ".metadata"``. If this is not possible,
        the variants metadata will be read from the BGEN file itself. If
        ``filepath + ".metadata"`` does not exist, it will try to create one
        with the same name to speed up reads. If ``False``, variants metadata
        will be read only from the BGEN file. If a file path is given instead,
        it assumes that the specified metadata file is valid and readable and
        therefore it will read variants metadata from that file only. Defaults
        to ``True``.
    sample_file : str, optional
        A sample file in `GEN format <http://www.stats.ox.ac.uk/~marchini/software/gwas/file_format.html>`_.
        If sample_file is provided, sample IDs are read from this file. Otherwise, it
        reads from the BGEN file itself if present. Defaults to ``None``.

    Returns
    -------
    dict
        variants : Variant position, chromossomes, RSIDs, etc.
        samples : Sample identifications.
        genotype : Array of genotype references.

    Notes
    -----
    Metadata files can speed up subsequent reads tremendously. But often the user does
    not have write permission for the default metadata file location
    ``filepath + ".metadata"``. We thus provide the
    :func:`bgen_reader.create_metadata_file` function for creating one at the
    given path.
    """
    import bgen_reader

    return bgen_reader.read_bgen(
        filepath,
        size=size,
        verbose=verbose,
        metadata_file=metadata_file,
        sample_file=sample_file,
    )


def convert_to_dosage(p, nalleles, ploidy):
    import bgen_reader

    return bgen_reader.convert_to_dosage(p, nalleles, ploidy)
