from __future__ import print_function

import sys
from colors import color


def oprint(*args, **kwargs):
    r"""Info print piped to stdout.

    Flush immediately.
    """
    print(*args, file=sys.stdout, flush=True, **kwargs)


def wprint(*args, **kwargs):
    r"""Warning print piped to stderr.

    Flush immediately.
    """
    msg = str(*args)
    print(color(msg, style='bold'), file=sys.stderr, flush=True, **kwargs)


def eprint(*args, **kwargs):
    r"""Error print piped to stderr.

    Flush immediately
    """
    msg = str(*args)
    print(
        color(msg, fg='red', style='bold'),
        file=sys.stderr,
        flush=True,
        **kwargs)
