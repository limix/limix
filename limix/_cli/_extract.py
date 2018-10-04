import click
import limix


@click.command()
@click.pass_context
@click.argument("filepath", type=click.Path(exists=True))
def extract(ctx, filepath):
    """Extract a file."""
    limix.sh.extract(filepath, verbose=ctx.obj["verbose"])
