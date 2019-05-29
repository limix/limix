import click

import limix
from ._misc import verbose_option


@click.command()
@click.argument("filepath", type=click.Path(exists=False))
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
@verbose_option
def see(filepath, show_chunks, header, verbose):
    """
    Show an overview of file or group of files.

    It supports HDF5, CSV,
    """
    from limix.io._fetch import parse_fetch_spec

    spec = parse_fetch_spec(filepath)
    filetype = spec.filetype
    filepath = spec.filepath

    if filetype == "unknown":
        print("Unknown file type or file path not reachable: `%s`." % filepath)

    # if filetype == "guess":
    #     filepath, filetype = limix.io.detect_file_type(filepath)

    if filetype == "hdf5":
        limix.io.hdf5._see(filepath, show_chunks=show_chunks)

    elif filetype == "csv":
        limix.io.csv._see(filepath, verbose=verbose, header=header)

    elif filetype == "plink2-rel":
        raise NotImplementedError
        # limix.io.plink._see_rel(filepath + ".rel", filepath + ".rel", verbose)
        # limix.plot.show()

    elif filetype == "grm.raw":
        limix.io.plink._see_kinship(filepath, verbose)
        limix.plot.show()

    elif filetype == "bed":
        limix.io.plink._see_bed(filepath, verbose)

    elif filetype == "bimbam-pheno":
        limix.io.bimbam._see_phenotype(filepath, verbose)

    elif filetype == "npy":
        limix.io.npy._see(filepath, verbose)

    elif filetype == "image":
        r = limix.plot.image(filepath)
        limix.plot.show()
        return r
