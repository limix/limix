import click
from ._misc import verbose_option


@click.command()
@click.argument("url")
@click.option(
    "--dest", help="Destination folder.", default=None, type=click.Path(exists=True)
)
@verbose_option
def download(url, dest, verbose):
    """
    Download file from the specified URL.
    """
    import limix

    limix.sh.download(url, dest=dest, verbose=verbose)
