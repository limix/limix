import click

import limix


@click.command()
@click.pass_context
@click.argument("url")
@click.option(
    "--dest", help="Destination path.", default=None, type=click.Path(exists=True)
)
@click.option(
    "--verbose/--quiet", "-v/-q", help="Enable or disable verbose mode.", default=True
)
def download(ctx, url, dest, verbose):
    """Download file from the specified URL."""
    limix.sh.download(url, dest=dest, verbose=verbose)
