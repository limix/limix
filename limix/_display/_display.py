import sys
from time import time
from ._core import pprint, bold, blue, red, width, wrap_text


def banner():
    from limix import __version__

    pyver = sys.version.split("\n")[0].strip()
    return "Running Limix {} using Python {}.".format(__version__, pyver)


def add_title_header(title, df):
    msg = repr(df)
    k = msg.find("\n") - len(title) - 2
    left = ("-" * (k // 2)) + " "
    right = " " + ("-" * (k // 2 + k % 2))
    out = left + title + right + "\n"
    out += msg
    return out


class session_line(object):
    r"""Print the elapsed time after the execution of a block of code."""

    def __init__(self, desc="Running... ", disable=False):
        self._disable = disable
        self._tstart = None
        self._desc = desc
        self.elapsed = None

    def __enter__(self):
        self._tstart = time()
        if not self._disable:
            sys.stdout.write(self._desc)
            sys.stdout.flush()
        return self

    def __exit__(self, *args, **_):
        from humanfriendly import format_timespan

        self.elapsed = time() - self._tstart
        if not self._disable:
            print("done (%s)." % format_timespan(self.elapsed))
            sys.stdout.flush()


class session_block(object):
    """Print session block: session start and session end."""

    def __init__(self, session_name, disable=False):
        self._session_name = session_name
        self._start = None
        self._disable = disable

    def __enter__(self):
        self._start = time()
        msg = " {} session starts ".format(self._session_name)
        if not self._disable:
            msg = wrap_text(msg, width())
            pprint(bold(blue(msg)))

    def __exit__(self, exception_type, exception_value, traceback):
        elapsed = time() - self._start
        fail = exception_type is not None

        if fail:
            msg = " {} session fails in {:.2f} seconds "
            color = red
        else:
            msg = " {} session ends in {:.2f} seconds "
            color = blue

        msg = msg.format(self._session_name, elapsed)
        if not self._disable:
            msg = wrap_text(msg, width())
            pprint(bold(color(msg)))
