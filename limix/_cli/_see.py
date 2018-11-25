import limix
import click


@click.command()
@click.pass_context
@click.argument("filepath")
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
@click.option(
    "--verbose/--quiet", "-v/-q", help="Enable or disable verbose mode.", default=True
)
def see(ctx, filepath, show_chunks, header, verbose):
    """ Show an overview of multiple file types. """
    from ._spec import parse_fetch_spec

    spec = parse_fetch_spec(filepath)
    filetype = spec["filetype"]
    filepath = spec["filepath"]

    if filetype == "unknown":
        print("Unknown file type or file path not reachable: `%s`." % filepath)

    if filetype == "guess":
        filepath, filetype = limix.io.detect_file_type(filepath)

    if filetype == "hdf5":
        limix.io.hdf5.see(filepath, show_chunks=show_chunks)

    elif filetype == "csv":
        limix.io.csv.see(filepath, verbose=verbose, header=header)

    elif filetype == "grm.raw":
        r = limix.io.plink.see_kinship(filepath, verbose)
        limix.plot.show()
        return r

    elif filetype == "bed":
        limix.io.plink.see_bed(filepath, verbose)

    elif filetype == "bimbam-pheno":
        limix.io.bimbam.see_phenotype(filepath, verbose)

    elif filetype == "npy":
        limix.io.npy.see(filepath, verbose)

    elif filetype == "image":
        r = limix.plot.image(filepath)
        limix.plot.show()
        return r
