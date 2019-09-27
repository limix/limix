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
    from .conftest import setup_tests_baseline

    setup_tests_baseline()

    args = [
        "--doctest-plus",
        "--doctest-plus-rtol=1e-02",
        "--doctest-plus-atol=1e-02",
        "--doctest-modules",
    ]
    if not verbose:
        args += ["--quiet"]

    args += ["--pyargs", __name__.split(".")[0]]

    return __import__("pytest").main(args)
