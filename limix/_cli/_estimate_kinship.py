import click
import limix


@click.command()
@click.pass_context
@click.argument("input-file", type=click.Path(exists=True))
@click.option(
    "--output-file",
    help="Specify the output file path.",
    default=None,
    type=click.Path(exists=False),
)
@click.option(
    "--filetype", help="Specify the file type instead of guessing.", default="guess"
)
def estimate_kinship(ctx, input_file, output_file, filetype):
    """Estimate a kinship matrix."""

    if filetype == "guess":
        filetype = limix.io.detect_file_type(input_file)

    if ctx.obj["verbose"]:
        print("Detected file type: {}".format(filetype))

    if filetype == "bgen":
        G = limix.io.bgen.fetch_dosage(input_file, verbose=ctx.obj["verbose"])
    elif filetype == "bed":
        G = limix.io.plink.fetch_dosage(input_file, verbose=ctx.obj["verbose"])
    else:
        print("Unknown file type: %s" % input_file)

    K = limix.stats.linear_kinship(G, verbose=ctx.obj["verbose"])

    if output_file is None:
        output_file = input_file + ".npy"

    output_file, oft = limix.io.detect_file_type(output_file)

    if oft == "npy":
        limix.io.npy.save_kinship(output_file, K, verbose=ctx.obj["verbose"])
    else:
        print("Unknown output file type: %s" % output_file)
