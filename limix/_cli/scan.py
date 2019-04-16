import click

from ._click import OrderedCommand


@click.command(cls=OrderedCommand)
@click.pass_context
@click.argument("trait")
@click.argument("genotype")
@click.option(
    "--covariate",
    help="Specify the file path to a file containing the covariates.",
    default=None,
)
@click.option(
    "--kinship",
    help="Specify the file path to a file containing the kinship matrix.",
    default=None,
)
@click.option(
    "--lik",
    help=(
        "Specify the type of likelihood that will described the residual error "
        "distribution. It can be 'normal', 'bernoulli', 'binomial', or 'poisson'."
    ),
    default="normal",
)
@click.option(
    "--drop-maf",
    help=(
        "Drop out candidates having a minor allele frequency below the provided "
        "threshold."
    ),
)
@click.option(
    "--drop-missing",
    help=("Drop out samples, candidates, or covariates with missing values."),
    multiple=True,
)
@click.option(
    "--where",
    help=(
        "Filtering expression to select which phenotype, genotype loci, and covariates"
        " to use in the analysis. The syntax is <TARGET>:<COND>,"
        " where <COND> is the first argument of the method `DataArray.where`."
    ),
    multiple=True,
)
@click.option(
    "--impute",
    help=(
        "Impute missing values for phenotype, genotype, and covariate. The syntax is"
        " <TARGET>:<DIM>:<METHOD>, where <METHOD> can be 'mean'. Defaults to 'mean'"
    ),
    multiple=True,
)
@click.option(
    "--normalize",
    help=(
        "Normalize phenotype, genotype, and covariate. The syntax is "
        "<TARGET>:<DIM>:<METHOD>, where <METHOD> can be 'gaussianize'"
        " or 'mean_std'. Defaults to 'gaussianize'."
    ),
    multiple=True,
)
@click.option(
    "--output-dir", help="Specify the output directory path.", default="output"
)
@click.option(
    "--verbose/--quiet", "-v/-q", help="Enable or disable verbose mode.", default=True
)
@click.option(
    "--dry-run/--no-dry-run",
    help="Perform a trial run with no scan taking place.",
    default=False,
)
def scan(
    ctx, trait, genotype, covariate, kinship, lik, output_dir, verbose, dry_run, **_
):
    """ Single-variant association testing via mixed models.

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
    import sys
    from os import makedirs
    from os.path import abspath, exists, join
    import traceback
    from limix._display import session_block, banner, session_line, indent, print_exc
    from limix.qtl import scan
    from limix.io import fetch
    from .pipeline import Pipeline
    from limix._data import conform_dataset
    from .preprocess import impute as impute_func
    from .preprocess import normalize as normalize_func
    from .preprocess import where as where_func
    from .preprocess import drop_missing, drop_maf

    print(banner())

    if ctx.obj is None:
        ctx.obj = {"preprocess": []}

    output_dir = abspath(output_dir)
    if not dry_run:
        if not exists(output_dir):
            makedirs(output_dir, exist_ok=True)

    def _print_data_array(arr, verbose):
        if verbose:
            print("\n{}\n".format(indent(_clean_data_array_repr(arr))))

    data = {"y": None, "G": None, "K": None}

    data["y"] = fetch("trait", trait, verbose)
    _print_data_array(data["y"], verbose)

    data["G"] = fetch("genotype", genotype, verbose)
    _print_data_array(data["G"], verbose)

    if covariate is not None:
        data["M"] = fetch("covariate", covariate, verbose)
        _print_data_array(data["M"], verbose)

    if kinship is not None:
        data["K"] = fetch("kinship", kinship, verbose)
        _print_data_array(data["K"], verbose)

    with session_line("Matching samples... "):
        data = conform_dataset(**data)
    data = {k: v for k, v in data.items() if v is not None}

    if data["y"].sample.size == 0:
        raise RuntimeError(
            "Exiting early because there is no sample left after matching samples."
            + " Please, check your sample ids."
        )

    oparams = _ordered_params(ctx)

    with session_block("preprocessing", disable=not verbose):
        pipeline = Pipeline(data)
        preproc_params = [
            i
            for i in oparams
            if i[0] in ["impute", "normalize", "where", "drop_missing", "drop_maf"]
        ]

        for p in preproc_params:
            if p[0] == "where":
                pipeline.append(where_func, "where", p[1])
            elif p[0] == "normalize":
                pipeline.append(normalize_func, "normalize", p[1])
            elif p[0] == "impute":
                pipeline.append(impute_func, "impute", p[1])
            elif p[0] == "drop_maf":
                pipeline.append(drop_maf, "drop-maf", p[1])
            elif p[0] == "drop_missing":
                pipeline.append(drop_missing, "drop-missing", p[1])

        data = pipeline.run()

    if dry_run:
        print("Exiting early because of dry-run option.")
        return

    if "K" not in data:
        data["K"] = None
    try:
        res = scan(
            data["G"], data["y"], lik=lik, K=data["K"], M=data["M"], verbose=verbose
        )
    except Exception as e:
        print_exc(traceback.format_stack(), e)
        sys.exit(1)

    with session_line("Saving results to `{}`... ".format(output_dir)):
        res.to_csv(join(output_dir, "null.csv"), join(output_dir, "alt.csv"))


def _clean_data_array_repr(arr):
    txt = str(arr)
    txt = txt.replace("xarray.DataArray ", "")
    txt = txt.replace("object ", "")
    txt = txt.replace("int64 ", "")
    txt = txt.replace("<U5 ", "")
    txt = txt.replace("dask.array", "array")
    return txt


def _ordered_params(ctx):
    args_seq = []
    for p in ctx.param_order:
        for opt, val in ctx.params.items():
            if p.name == opt:
                if isinstance(val, tuple):
                    v = val[0]
                    val = val[1:]
                else:
                    v = val
                args_seq.append((opt, v))
                if isinstance(val, tuple):
                    if len(val) == 0:
                        del ctx.params[opt]
                    else:
                        ctx.params[opt] = v
                else:
                    del ctx.params[opt]
                break
            pass
    return args_seq
