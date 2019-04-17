def test(verbose=True):
    """
    Run tests to verify this package's integrity.

    Parameters
    ----------
    verbose : bool
        ``True`` to show diagnostic. Defaults to ``True``.

    Returns
    -------
    int
        Exit code: ``0`` for success.
    """
    from .conftest import pandas_format

    pandas_format()

    args = [
        "--doctest-plus",
        "--doctest-plus-rtol=1e-04",
        "--doctest-plus-atol=1e-04",
        "--doctest-modules",
    ]
    if not verbose:
        args += ["--quiet"]

    args += ["--pyargs", __name__.split(".")[0]]

    return __import__("pytest").main(args)
