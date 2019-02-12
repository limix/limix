import click


@click.command()
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
        " to use in the analysis. The syntax is `<TARGET>: <COND>`,"
        " where <COND> is the first argument of the method `DataArray.where`."
    ),
    multiple=True,
)
@click.option(
    "--filter-missing",
    help=("Drop out samples, candidates, or covariates with missing values."),
    multiple=True,
)
@click.option(
    "--filter-maf",
    help=(
        "Drop out candidates having a minor allele frequency below the provided "
        "threshold."
    ),
)
@click.option(
    "--impute",
    help=("Impute missing values for phenotype, genotype, and covariate."),
    multiple=True,
)
@click.option(
    "--normalize", help=("Normalize phenotype, genotype, and covariate."), multiple=True
)
@click.option(
    "--output-dir", help="Specify the output directory path.", default="output"
)
@click.option(
    "--verbose/--quiet", "-v/-q", help="Enable or disable verbose mode.", default=True
)
@click.option(
    "--model",
    help=("Specify the statistical model to perform the scan."),
    default="single-trait",
    type=click.Choice(["single-trait", "struct-lmm"]),
)
@click.option(
    "--dry-run/--no-dry-run",
    help="Perform a trial run with no scan taking place.",
    default=False,
)
def scan(
    ctx,
    phenotypes_file,
    genotype_file,
    covariates_file,
    kinship_file,
    lik,
    filter,
    filter_missing,
    filter_maf,
    impute,
    normalize,
    output_dir,
    verbose,
    model,
    dry_run,
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
    from os import makedirs
    from os.path import abspath, exists
    from limix._display import session_block, banner, session_line, indent
    from .._fetch import fetch
    from .pipeline import Pipeline
    from limix._data import conform_dataset
    from .._preprocess import process_filter, process_normalize

    print(banner())

    output_dir = abspath(output_dir)
    if not exists(output_dir):
        makedirs(output_dir, exist_ok=True)

    def _print_data_array(arr, verbose):
        if verbose:
            print("\n{}\n".format(indent(_clean_data_array_repr(arr))))

    data = {"y": None, "G": None, "K": None}

    data["y"] = fetch(phenotypes_file, "trait", verbose)
    _print_data_array(data["y"], verbose)

    data["G"] = fetch(genotype_file, "genotype", verbose)
    _print_data_array(data["G"], verbose)

    if kinship_file is not None:
        data["K"] = fetch(kinship_file, "kinship", verbose)
        _print_data_array(data["K"], verbose)

    with session_line("Matching samples... "):
        data = conform_dataset(**data)
    data = {k: v for k, v in data.items() if v is not None}

    if data["y"].sample.size == 0:
        raise RuntimeError(
            "Exiting early because there is no sample left after matching samples."
            + " Please, check your sample ids."
        )

    with session_block("preprocessing", disable=not verbose):
        pipeline = Pipeline(data)

        for spec in filter:
            for target in data.keys():
                pipeline.append(process_filter, spec, target)

        for spec in normalize:
            for target in data.keys():
                pipeline.append(process_normalize, spec, target)

        pipeline.run()

    # with session_block("preprocessing", disable=not verbose):
    #     data = _preprocessing(
    #         data, filter, filter_missing, filter_maf, impute, normalize, verbose
    #     )

    # if dry_run:
    #     print("Exiting early because of dry-run option.")
    #     return

    # if "K" not in data:
    #     data["K"] = None
    # try:
    #     model = limix.qtl.st_scan(
    #         data["G"], data["y"], lik, K=data["K"], verbose=verbose
    #     )
    # except Exception as e:
    #     print_exc(traceback.format_stack(), e)
    #     sys.exit(1)

    # with session_line("Saving results to `{}`... ".format(output_dir)):
    #     model.to_csv(join(output_dir, "null.csv"), join(output_dir, "alt.csv"))


def _clean_data_array_repr(arr):
    txt = str(arr)
    txt = txt.replace("xarray.DataArray ", "")
    txt = txt.replace("object ", "")
    txt = txt.replace("int64 ", "")
    txt = txt.replace("<U5 ", "")
    txt = txt.replace("dask.array", "array")
    return txt
