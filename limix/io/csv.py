def read(filename, sep=None, header=True, verbose=True):
    """
    Read a CSV file.

    Parameters
    ----------
    filename : str
        Path to a CSV file.
    sep : str
        Separator. ``None`` triggers auto-detection. Defaults to ``None``.
    header : bool
        ``True`` for file with a header; ``False`` otherwise. Defaults
        to ``True``.
    verbose : bool
        `True` for progress information; `False` otherwise.

    Returns
    -------
    data : dask dataframes

    Examples
    --------
    .. doctest::

        >>> from limix.io.csv import read
        >>> from limix import file_example
        >>>
        >>> with file_example("data.csv") as filepath:
        ...     df = read(filepath, verbose=False)
        ...     print(df)  # doctest: +FLOAT_CMP
           pheno   attr1 attr2 attr3
        0    sex  string    10     a
        1   size   float    -3     b
        2  force     int     f     c
    """
    from dask.dataframe import read_csv as dask_read_csv
    from pandas import read_csv as pandas_read_csv
    from .._display import session_line

    if sep is None:
        sep = _infer_separator(filename)

    header = 0 if header else None

    with session_line("Reading {}... ".format(filename), disable=not verbose):

        if _is_large_file(filename):
            df = dask_read_csv(filename, sep=sep, header=header)
        else:
            df = pandas_read_csv(filename, sep=sep, header=header)

    if len(df.columns) > 0:
        if df.columns[0] == "Unnamed: 0":
            df = df.set_index("Unnamed: 0")
            df.index.name = None

    return df


def _see(filepath, header, verbose=True):
    """
    Shows a human-friendly representation of a CSV file.

    Parameters
    ----------
    filepath : str
        CSV file path.
    header : bool
        ``True`` for parsing the header; ``False`` otherwise.
    verbose : bool
        ``True`` for verbose; ``False`` otherwise.

    Returns
    -------
    str
        CSV representation.
    """
    from pandas import read_csv
    from .._display import session_line

    if header:
        header = 0
    else:
        header = None

    with session_line(desc="Reading %s... " % filepath, disable=not verbose):
        sep = _infer_separator(filepath)
        msg = read_csv(filepath, sep=sep, header=header).head()

    print(msg)


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

    with open(fn, "r") as f:
        line = _remove_repeat(f.readline())
        counter = _count(set(line), line)

        for _ in range(nmax - 1):
            line = _remove_repeat(f.readline())
            if len(line) == 0:
                break
            c = _count(set(counter.keys()), line)
            _update(counter, c)
            if len(counter) == 1:
                return next(iter(counter.keys()))

    for c in set([",", "\t", " "]):
        if c in counter:
            return c

    counter = list(counter.items())
    if len(counter) == 0:
        return None

    counter = sorted(counter, key=lambda kv: kv[1])
    return counter[-1][0]


def _remove_repeat(s):
    from re import sub

    return sub(r"(.)\1+", r"\1", s)


def _is_large_file(filepath):
    import os

    large = 1024 * 1024 * 100
    return os.path.getsize(filepath) >= large
