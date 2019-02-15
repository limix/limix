import click
from ._click import OrderedCommand


@click.command(cls=OrderedCommand)
@click.pass_context
def st_sscan():
    """ Struct-LMM. """
    pass
