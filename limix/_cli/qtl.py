import click

from .qtl_scan import scan


# @click.command(cls=OrderedCommand)
@click.group()
# @click.option(
#     "--model",
#     help=("Specify the statistical model to perform the scan."),
#     default="single-trait",
#     type=click.Choice(["single-trait", "struct-lmm"]),
# )
# @click.option(
#     "--lik",
#     help=(
#         "Specify the type of likelihood that will described"
#         " the residual error distribution."
#     ),
#     default="normal",
# )
def qtl():
    """
    Perform genome-wide association scan.
    """
    pass


qtl.add_command(scan)
