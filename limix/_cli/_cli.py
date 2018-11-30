import limix
import click
from ._see import see
from ._estimate_kinship import estimate_kinship
from ._download import download
from ._extract import extract
from ._scan import scan
from ._remove import remove


@click.group(name="limix", context_settings=dict(help_option_names=["-h", "--help"]))
@click.pass_context
@click.version_option(version=limix.__version__)
def cli(ctx):
    pass


cli.add_command(see)
cli.add_command(estimate_kinship)
cli.add_command(download)
cli.add_command(extract)
cli.add_command(scan)
cli.add_command(remove)
