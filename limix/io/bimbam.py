def read_phenotype(filepath):
    r"""Read a BIMBAM phenotype file.

    Parameters
    ----------
    filepath : str
        File path.

    Returns
    -------
    :class:`pandas.DataFrame`
        DataFrame representation of the file.

    Examples
    --------
    .. doctest::

        >>> import limix
        >>> from limix.example import file_example
        >>>
        >>> with file_example("ex0/phenotype.gemma") as filepath:
        ...     print(limix.io.bimbam.read_phenotype(filepath))
                 0        1        2
        0  1.20000 -0.30000 -1.50000
        1      nan  1.50000  0.30000
        2  2.70000  1.10000      nan
        3 -0.20000 -0.70000  0.80000
        4  3.30000  2.40000  2.10000
    """
    from pandas import read_csv

    df = read_csv(filepath, sep=r"\s+", header=None)

    return df


def see_phenotype(filepath):
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
    from pandas import read_csv

    df = read_csv(filepath, sep=r"\s+", header=None)
    print(df.head().to_string(show_dimensions=True))
