import sys
from limix.display import bold, red, pprint


def print_exc(stack, e):

    print("Traceback (most recent call last):")
    sys.stdout.write("".join(stack))
    errmsg = "{}: {}".format(e.__class__.__name__, str(e))
    pprint(bold(red(errmsg)))
