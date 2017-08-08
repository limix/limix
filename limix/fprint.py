from __future__ import print_function

import sys


def oprint(*args, **kwargs):
    print(*args, file=sys.stdout, flush=True, **kwargs)


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, flush=True, **kwargs)
