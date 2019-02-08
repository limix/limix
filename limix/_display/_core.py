r"""Retrieve information about terminal/frontend.

Acknowledgment
--------------
- Pandas Python package for providing code for such a functionality.
"""
import re
import warnings

from .._config import config

_tags = ["bold", "green", "blue", "red"]
_color = {"blue": "#0C68C7", "green": "#19CB00", "red": "#FF3534"}


def pprint(txt):
    """ Print rich text. """
    try:
        from IPython.display import display

        display(_RichText(txt))
    except Exception:
        print(_RichText(txt))


def add_title_header(title, df):
    msg = repr(df)
    k = msg.find("\n") - len(title) - 2
    left = ("-" * (k // 2)) + " "
    right = " " + ("-" * (k // 2 + k % 2))
    out = left + title + right + "\n"
    out += msg
    return out


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


def red(txt):
    """Red color font."""
    return "[red]" + txt + "[/red]"


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


def wrap_text(msg, width, sym="="):
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
    try:
        from blessings import Terminal
    except ImportError:
        return None

    try:
        term = Terminal()
    except ValueError:
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

        return "<pre>{}</pre>".format(txt)


def _terminal_format(txt):
    term = _get_terminal()

    for tag in _tags:
        r = _compile(tag)
        if term is None:
            by = r"\1"
        else:
            by = getattr(term, tag)(r"\1")
        txt = r.sub(by, txt)

    return txt


def _plain_format(txt):

    for tag in _tags:
        r = _compile(tag)
        txt = r.sub("\\1", txt)

    return txt
