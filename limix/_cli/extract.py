import click

import limix


@click.command()
@click.pass_context
@click.argument("filepath", type=click.Path(exists=True))
@click.option(
    "--verbose/--quiet", "-v/-q", help="Enable or disable verbose mode.", default=True
)
def extract(ctx, filepath, verbose):
    """Extract a file."""
    limix.sh.extract(filepath, verbose=verbose)
