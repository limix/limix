import click

from ._misc import get_version
from .download import download
from .estimate_kinship import estimate_kinship
from .extract import extract
from .qtl import qtl
from .remove import remove
from .see import see


@click.group(name="limix", context_settings=dict(help_option_names=["-h", "--help"]))
@click.pass_context
@click.version_option(get_version())
def cli(ctx):
    pass


cli.add_command(see)
cli.add_command(estimate_kinship)
cli.add_command(download)
cli.add_command(extract)
cli.add_command(qtl)
cli.add_command(remove)
