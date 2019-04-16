import sys

from ._display import bold, pprint, red


def print_exc(stack, e):
    """ Print exception with red and bold fonts. """

    print("Traceback (most recent call last):")
    sys.stdout.write("".join(stack))
    errmsg = "{}: {}".format(e.__class__.__name__, str(e))
    pprint(bold(red(errmsg)))
