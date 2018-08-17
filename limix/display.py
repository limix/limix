r"""Retrieve information about terminal/frontend.

Acknowledgment
--------------
- Pandas Python package for providing code for such a functionality.
"""
import io
import re
import sys
import warnings
from contextlib import redirect_stdout
from time import time

from ._config import config

_tags = ["bold", "green", "blue"]
_color = {"blue": "#0C68C7", "green": "#19CB00"}


def pprint(txt):
    """Print rich text."""
    try:
        from IPython.display import display

        display(_RichText(txt))
    except Exception:
        print(_RichText(txt))


def format_richtext(txt):
    """Format rich text."""
    return _RichText(txt)


def display(objs):
    """Frontend-aware object display."""
    try:
        from IPython.display import display

        display(objs)
    except Exception:
        print(objs)


def bold(txt):
    """Bold font."""
    return "[bold]" + txt + "[/bold]"


def green(txt):
    """Green color font."""
    return "[green]" + txt + "[/green]"


def blue(txt):
    """Blue color font."""
    return "[blue]" + txt + "[/blue]"


def width():
    """Display number of columns."""
    if _in_interactive_session():
        if _in_ipython_frontend():
            return config["display.fallback_width"]
    try:
        term = _get_terminal()
        if term is not None:
            return term.width
    except ValueError as e:
        warnings.warn(e)

    return config["display.fallback_width"]


class timer_text(object):
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

    def __exit__(self, *args, **kwargs):
        from humanfriendly import format_timespan

        self.elapsed = time() - self._tstart
        if not self._disable:
            print("done (%s)." % format_timespan(self.elapsed))
            sys.stdout.flush()


class session_text(object):
    """Print session block: session start and session end."""

    def __init__(self, session_name, disable=False):
        self._session_name = session_name
        self._start = None
        self._disable = disable

    def __enter__(self):
        self._start = time()
        msg = " {} session starts ".format(self._session_name)
        if not self._disable:
            msg = _msg_wrap(msg, width())
            pprint(bold(blue(msg)))

    def __exit__(self, *_):
        elapsed = time() - self._start
        msg = " {} session ends in {:.2f} seconds "
        msg = msg.format(self._session_name, elapsed)
        if not self._disable:
            msg = _msg_wrap(msg, width())
            pprint(bold(blue(msg)))


def _msg_wrap(msg, width, sym="="):
    if width is None:
        width = 79
    width -= len(msg)
    if width <= 0:
        return msg
    left = width // 2
    right = width - left
    return (sym * left) + msg + (sym * right)


def _in_ipython_frontend():
    """Check if we're inside an an IPython zmq frontend."""
    try:
        ip = get_ipython()  # noqa
        return "zmq" in str(type(ip)).lower()
    except Exception:
        pass
    return False


def _is_terminal():
    """Detect if Python is running in a terminal.

    Returns
    -------
    bool
        ``True`` if Python is running in a terminal; ``False`` otherwise.
    """
    try:
        ip = get_ipython()
    except NameError:  # assume standard Python interpreter in a terminal
        return True
    else:
        if hasattr(ip, "kernel"):  # IPython as a Jupyter kernel
            return False
        else:  # IPython in a terminal
            return True


def _in_interactive_session():
    """Check if we're running in an interactive shell.

    Returns
    -------
    ``True`` if running under python/ipython interactive shell; ``False`` otherwise.
    """

    def check_main():
        import __main__ as main

        return not hasattr(main, "__file__")

    try:
        return __IPYTHON__ or check_main()  # noqa
    except Exception:
        return check_main()


def _get_terminal():
    from blessings import Terminal

    try:
        term = Terminal()
    except ValueError as e:
        return None
    return term


def _compile(tag):
    expr = r"\[{tag}\](.*?)\[\/{tag}\]".format(tag=tag)
    return re.compile(expr, re.MULTILINE | re.DOTALL)


class _RichText(object):
    def __init__(self, text):
        self._text = text

    def __repr__(self):
        if _is_terminal() or _in_ipython_frontend():
            return _terminal_format(self._text)
        return _plain_format(self._text)

    def _repr_html_(self):
        txt = self._text

        for tag in _tags:
            r = _compile(tag)
            if tag == "bold":
                txt = r.sub("<b>\\1</b>".format(tag), txt)
            else:
                expr = "<span style='color:{}'>\\1</span>".format(_color[tag])
                txt = r.sub(expr, txt)

        # return txt.replace("\n", "\n<br>")
        return "<pre>{}</pre>".format(txt)


def _terminal_format(txt):
    from limix.display import _get_terminal

    term = _get_terminal()

    for tag in _tags:
        r = _compile(tag)
        txt = r.sub(getattr(term, tag)("\\1"), txt)

    return txt


def _plain_format(txt):

    for tag in _tags:
        r = _compile(tag)
        txt = r.sub("\\1", txt)

    return txt
