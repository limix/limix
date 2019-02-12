import click
import limix


@click.command()
@click.pass_context
@click.argument("filepath", type=click.Path(exists=True))
@click.option(
    "--verbose/--quiet", "-v/-q", help="Enable or disable verbose mode.", default=True
)
@click.option(
    "--dest", help="Destination path.", default=".", type=click.Path(exists=True)
)
def extract(ctx, filepath, dest, verbose):
    """Extract a file."""
    limix.sh.extract(filepath, dest=dest, verbose=verbose)
