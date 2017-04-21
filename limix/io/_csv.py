from dask.dataframe import read_csv as _read_csv


def read_csv(filename):
    r"""Read a CSV file.

    Args:
        filename (str): path to a CSV file.

    Returns:
        dask.DataFrame: a dataframe representation.

    Example
    -------

        .. doctest::

            >>> from limix.io import read_csv
            >>> from limix.io.examples import csv_file_example
            >>>
            >>> df = read_csv(csv_file_example())
            >>> print(df.compute()) #doctest: +NORMALIZE_WHITESPACE
               pheno   attr1 attr2 attr3
            0    sex  string    10     a
            1   size   float    -3     b
            2  force     int     f     c

    """
    df = _read_csv(filename)
    df.set_index(df.columns[0], inplace=True)
    return df
