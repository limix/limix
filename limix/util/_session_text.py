from time import time
from blessings import Terminal
import warnings


class MockTerminal(object):
    @property
    def width(self):
        return 88

    def _return_msg(self, msg):
        return msg

    def __getattr__(self, name):
        return self._return_msg


class session_text(object):
    def __init__(self, session_name, disable=False):
        self._session_name = session_name
        self._start = None
        try:
            self._term = Terminal()
        except ValueError as e:
            warnings.warn(str(e))
            self._term = MockTerminal()
        self._disable = disable

    def __enter__(self):
        self._start = time()
        msg = " {} session starts ".format(self._session_name)
        if not self._disable:
            msg = _msg_wrap(msg, self._term.width)
            print(self._term.bold(self._term.bright_blue(msg)))

    def __exit__(self, *_):
        elapsed = time() - self._start
        msg = " {} session ends in {:.2f} seconds "
        msg = msg.format(self._session_name, elapsed)
        if not self._disable:
            msg = _msg_wrap(msg, self._term.width)
            print(self._term.bold(self._term.bright_blue(msg)))


def _msg_wrap(msg, width, sym="="):
    if width is None:
        width = 79
    width -= len(msg)
    if width <= 0:
        return msg
    left = width // 2
    right = width - left
    return (sym * left) + msg + (sym * right)
