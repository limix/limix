from os import chdir, getcwd
from os.path import abspath, basename, dirname, realpath


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

    pkgname = basename(dirname(realpath(__file__)))

    p = __import__(pkgname).__path__[0]
    src_path = abspath(p)
    old_path = getcwd()
    chdir(src_path)

    try:
        return_code = __import__('pytest').main(['-q', '--doctest-modules'])
    finally:
        chdir(old_path)

    if return_code == 0:
        print("Congratulations. All tests have passed!")

    return return_code
