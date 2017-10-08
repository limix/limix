from __future__ import print_function

import sys


# TODO: add wprint, document it, and make sure everywhere is using those
def oprint(*args, **kwargs):
    print(*args, file=sys.stdout, flush=True, **kwargs)


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, flush=True, **kwargs)
