import limix
import click


@click.command()
@click.pass_context
@click.argument("filepath")
@click.option(
    "--filetype", help="Specify the file type instead of guessing.", default="guess"
)
@click.option(
    "--show_chunks",
    help="Chunks if datasets will be displayed, if available.",
    default="guess",
)
@click.option(
    "--header / --no-header",
    help="Parse header from CSV file. Defaults to false.",
    default=False,
)
def see(ctx, filepath, filetype, show_chunks, header):
    """Show an overview of multiple file types."""
    if filetype == "guess":
        filepath, filetype = limix.io.detect_file_type(filepath)

    if filetype == "hdf5":
        limix.io.hdf5.see(filepath, show_chunks=show_chunks)

    elif filetype == "csv":
        limix.io.csv.see(filepath, verbose=ctx.obj["verbose"], header=header)

    elif filetype == "grm.raw":
        r = limix.io.plink.see_kinship(filepath, verbose=ctx.obj["verbose"])
        limix.plot.show()
        return r

    elif filetype == "bed":
        limix.io.plink.see_bed(filepath, verbose=ctx.obj["verbose"])

    elif filetype == "bimbam-pheno":
        limix.io.bimbam.see_phenotype(filepath)

    elif filetype == "image":
        r = limix.plot.image(filepath)
        limix.plot.show()
        return r
    else:
        print("Unknown file type: %s" % filepath)
