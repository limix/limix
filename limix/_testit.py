def test(verbose=True):
    r"""Run tests to verify this package's integrity.

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
        "--doctest-modules",
        "--doctest-plus",
        "--doctest-plus-rtol=1e-03",
        "--doctest-plus-atol=1e-03",
    ]
    if not verbose:
        args += ["--quiet"]

    args += ["--pyargs", __name__.split(".")[0]]

    return __import__("pytest").main(args)
