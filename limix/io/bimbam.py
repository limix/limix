def read_phenotype(filepath, verbose=True):
    r"""Read a BIMBAM phenotype file.

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
        >>> from limix.example import file_example
        >>>
        >>> with file_example("ex0/phenotype.gemma") as filepath:
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
    from ..display import timer_text

    with timer_text("Reading `{}`... ".format(filepath), disable=not verbose):
        df = read_csv(filepath, sep=r"\s+", header=None)

    df.index = range(df.shape[0])
    df.index.name = "sample"
    df.columns = range(df.shape[1])
    df.columns.name = "trait"

    return df


def see_phenotype(filepath, verbose=True):
    """Shows a summary of a BIMBAM phenotype file.

    Parameters
    ----------
    filepath : str
        File path.

    Returns
    -------
    str
        File representation.
    """
    from ..display import dataframe_repr

    df = read_phenotype(filepath, verbose)

    print(dataframe_repr("Phenotypes", df))
