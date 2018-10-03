import click
from limix.io._detect import detect_file_type
from limix.io import fetch_phenotype, fetch_genotype


def _get_version():
    import limix

    return limix.__version__


def _show_plot():
    return
    import matplotlib as mpl

    if mpl.get_backend().lower() == "agg":
        return

    from limix import plot

    plt = plot.get_pyplot()
    plt.show()


@click.group(name="limix", context_settings=dict(help_option_names=["-h", "--help"]))
@click.pass_context
@click.option(
    "--verbose/--quiet", "-v/-q", help="Enable or disable verbose mode.", default=True
)
@click.version_option(version=_get_version())
def cli(ctx, verbose):
    ctx.obj = {}
    ctx.obj["verbose"] = verbose


@cli.command()
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
    from limix import io, plot

    if filetype == "guess":
        filepath, filetype = detect_file_type(filepath)

    if filetype == "hdf5":
        io.hdf5.see(filepath, show_chunks=show_chunks)

    elif filetype == "csv":
        io.csv.see(filepath, verbose=ctx.obj["verbose"], header=header)

    elif filetype == "grm.raw":
        r = io.plink.see_kinship(filepath, verbose=ctx.obj["verbose"])
        _show_plot()
        return r

    elif filetype == "bed":
        io.plink.see_bed(filepath, verbose=ctx.obj["verbose"])

    elif filetype == "bimbam-pheno":
        io.bimbam.see_phenotype(filepath)

    elif filetype == "image":
        r = plot.image(filepath)
        _show_plot()
        return r
    else:
        print("Unknown file type: %s" % filepath)


@cli.command()
@click.pass_context
@click.argument("url")
@click.option(
    "--dest", help="Destination path.", default=None, type=click.Path(exists=True)
)
@click.option("--force", help="Overwrite existing file if necessary.", is_flag=True)
def download(ctx, url, dest, force):
    """Download file from the specified URL."""
    import limix

    limix.sh.download(url, dest=dest, verbose=ctx.obj["verbose"], force=force)


@cli.command()
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
    from limix import io, stats

    if filetype == "guess":
        filetype = detect_file_type(input_file)

    if ctx.obj["verbose"]:
        print("Detected file type: {}".format(filetype))

    if filetype == "bgen":
        G = io.bgen.fetch_dosage(input_file, verbose=ctx.obj["verbose"])
    elif filetype == "bed":
        G = io.plink.fetch_dosage(input_file, verbose=ctx.obj["verbose"])
    else:
        print("Unknown file type: %s" % input_file)

    K = stats.linear_kinship(G, verbose=ctx.obj["verbose"])

    if output_file is None:
        output_file = input_file + ".npy"

    output_file, oft = detect_file_type(output_file)

    if oft == "npy":
        io.npy.save_kinship(output_file, K, verbose=ctx.obj["verbose"])
    else:
        print("Unknown output file type: %s" % output_file)


@cli.command()
@click.pass_context
@click.argument("filepath", type=click.Path(exists=True))
def extract(ctx, filepath):
    """Extract a file."""
    import limix

    limix.sh.extract(filepath, verbose=ctx.obj["verbose"])


@cli.command()
@click.pass_context
@click.argument("filepath", type=click.Path(exists=True))
def remove(ctx, filepath):
    """Remove a file."""
    import limix

    limix.sh.remove(filepath)


def _process_phenotype_filter(df, flt):
    import re

    flt = flt.strip()
    r = flt.replace("phenotype", "", 1).strip()
    m = re.match(r"^(\[.*\])+", r)
    if m is not None:
        slic = m.groups(0)[0]
        df = eval("df{}".format(slic))
    return df


_dispath_process_filter = {"phenotype": _process_phenotype_filter}


@cli.command()
@click.pass_context
@click.argument("phenotypes-file")
@click.argument("genotype-file")
@click.option(
    "--covariates-file",
    help="Specify the file path to a file containing the covariates.",
    default=None,
)
@click.option(
    "--kinship-file",
    help="Specify the file path to a file containing the kinship matrix.",
    default=None,
)
@click.option(
    "--lik",
    help=(
        "Specify the type of likelihood that will described"
        " the residual error distribution."
    ),
    default="normal",
)
@click.option(
    "--filter",
    help=(
        "Filtering expression to select which phenotype, genotype loci, and covariates"
        " to use in the analysis."
    ),
    multiple=True,
)
@click.option("--output-dir", help="Specify the output directory path.", default=None)
def scan(
    ctx,
    phenotypes_file,
    genotype_file,
    covariates_file,
    kinship_file,
    lik,
    filter,
    output_dir,
):
    """Perform genome-wide association scan.

    This analysis requires minimally the specification of one phenotype
    (PHENOTYPES_FILE) and genotype data (GENOTYPE_FILE).

    The --filter option allows for selecting a subset of the original dataset for
    the analysis. For example,

        --filter="genotype: (chrom == '3') & (pos > 100) & (pos < 200)"

    states that only loci of chromosome 3 having a position inside the range (100, 200)
    will be considered. The --filter option can be used multiple times in the same
    call. In general, --filter accepts a string of the form

        <DATA-TYPE>: <BOOL-EXPR>

    where <DATA-TYPE> can be phenotype, genotype, or covariate. <BOOL-EXPR> is a boolean
    expression involving row or column names. Please, consult `pandas.DataFrame.query`
    function from Pandas package for further information.
    \f

    Examples
    --------

    ... doctest::

        # First we perform a quick file inspection. This step is optional but is very
        # useful to check whether `limix` is able to read them and print out their
        # metadata.
        limix show phenotypes.csv
        limix show genotype.bgen
        limix show kinship.raw

        # We now perform the analysis, specifying the genotype loci and the phenotype
        # of interest.
        limix phenotypes.csv genotype.bgen --kinship-file=kinship.raw \
            --output-dir=results \
            --filter="phenotype: col == 'height'" \
            --filter="genotype: (chrom == '3') & (pos > 100) & (pos < 200)"
    """
    import limix

    pheno_filepath, pheno_type = detect_file_type(phenotypes_file)

    geno_filepath, geno_type = detect_file_type(genotype_file)

    y = fetch_phenotype(pheno_filepath, pheno_type)
    G = fetch_genotype(geno_filepath, geno_type)

    data = {"phenotype": y, "genotype": G}

    import pdb

    for filt in filter:
        target = _get_filter_target(filt)
        data[target] = _dispath_process_filter[target](data[target], filt)

    pdb.set_trace()

    r = limix.qtl.scan(data["genotype"], data["phenotype"], lik)
    # limix scan ex0/phenotype.gemma:bimbam-pheno genotype:be
    # scan(G, y, lik, K=None, M=None, verbose=True)


def _get_filter_target(flt):
    flt = flt.strip()
    a, b = flt.find("["), flt.find(":")
    s = _sign(a) * _sign(b)
    i = s * min(s * a, s * b)
    return flt[: i + 1][:-1].strip().lower()


def _sign(v):
    if v > 0:
        return 1
    if v < 0:
        return -1
    return 0
