import limix
import click


@click.command()
@click.pass_context
@click.argument("filepath", type=click.Path(exists=True))
def remove(ctx, filepath):
    """Remove a file."""

    limix.sh.remove(filepath)
