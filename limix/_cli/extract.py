import click

import limix
from ._misc import verbose_option


@click.command()
@click.argument("filepath", type=click.Path(exists=True))
@verbose_option
def extract(filepath, verbose):
    """
    Extract a file.
    """
    limix.sh.extract(filepath, verbose=verbose)
