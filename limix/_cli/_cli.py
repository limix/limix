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
@click.option(
    "--verbose/--quiet", "-v/-q", help="Enable or disable verbose mode.", default=True
)
@click.version_option(version=limix.__version__)
def cli(ctx, verbose):
    ctx.obj = {}
    ctx.obj["verbose"] = verbose


cli.add_command(see)
cli.add_command(estimate_kinship)
cli.add_command(download)
cli.add_command(extract)
cli.add_command(scan)
cli.add_command(remove)
