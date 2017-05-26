def read_csv(filename, sep=None, header=True):
    r"""Read a CSV file.

    Parameters
    ----------
    filename : str
        Path to a CSV file.
    sep : str
        Separator.

    Returns
    -------
    data : dask dataframes

    Examples
    --------
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
    from dask.dataframe import read_csv as _read_csv

    if sep is None:
        sep = _infer_separator(filename)

    header = 0 if header else None
    return _read_csv(filename, sep=sep, header=header)


def see(filepath):
    """Shows a human-friendly representation of a CSV file.

    Parameters
    ----------
    filepath : str
        CSV file path.

    Returns
    -------
    str
        CSV representation.
    """
    from pandas import read_csv
    sep = _infer_separator(filepath)
    print(read_csv(filepath, sep=sep).describe())


def _count(candidates, line):
    counter = {c: 0 for c in candidates}
    for i in line:
        if i in candidates:
            counter[i] += 1
    return counter


def _update(counter, c):
    for (k, v) in c.items():
        if counter[k] != v:
            del counter[k]


def _infer_separator(fn):
    nmax = 9

    with open(fn, 'r') as f:
        line = f.readline().strip()
        counter = _count(set(line), line)

        for i in range(nmax - 1):
            line = f.readline().strip()
            if len(line) == 0:
                break
            c = _count(set(counter.keys()), line)
            _update(counter, c)
            if len(counter) == 1:
                return counter.keys()[0]

    for c in set([',', '\t', ' ']):
        if c in counter:
            return c

    counter = list(counter.items())
    counter = sorted(counter, key=lambda kv: kv[1])
    return counter[-1][0]
