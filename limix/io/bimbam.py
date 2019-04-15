def read_phenotype(filepath, verbose=True):
    """
    Read a BIMBAM phenotype file.

    Parameters
    ----------
    filepath : str
        File path.

    Returns
    -------
    :class:`pandas.DataFrame`
        DataFrame representation of the file.
    verbose : bool
        `True` for progress information; `False` otherwise.

    Examples
    --------
    .. doctest::

        >>> import limix
        >>> from limix import file_example
        >>>
        >>> with file_example("phenotype.gemma") as filepath:
        ...     print(limix.io.bimbam.read_phenotype(filepath, verbose=False))
        trait         0        1        2
        sample
        0        1.20000 -0.30000 -1.50000
        1            nan  1.50000  0.30000
        2        2.70000  1.10000      nan
        3       -0.20000 -0.70000  0.80000
        4        3.30000  2.40000  2.10000

    Notes
    -----
    BIMBAM phenotype files do not explicitly define sample ids (nor trait ids) but their
    order of appearance is used to associate samples from different files. Therefore,
    we denote the first sample found in this file as ``0``, the second as ``1``, and so
    on. We apply the same reasoning for trait naming.
    """
    from pandas import read_csv
    from .._display import session_line

    with session_line("Reading `{}`... ".format(filepath), disable=not verbose):
        df = read_csv(filepath, sep=r"\s+", header=None)

    df.index.name = "sample"
    df.columns.name = "trait"

    return df


def _see_phenotype(filepath, verbose=True):
    """
    Shows a summary of a BIMBAM phenotype file.

    Parameters
    ----------
    filepath : str
        File path.

    Returns
    -------
    str
        File representation.
    """
    from .._display import add_title_header

    df = read_phenotype(filepath, verbose)

    print(add_title_header("Phenotypes", df))
