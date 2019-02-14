import click
from ._click import OrderedCommand


@click.command(cls=OrderedCommand)
@click.pass_context
@click.argument("trait")
@click.argument("genotype")
@click.argument("env")
@click.option(
    "--covariates",
    help="Specify the file path to a file containing the covariates.",
    default=None,
)
@click.option(
    "--kinship",
    help="Specify the file path to a file containing the kinship matrix.",
    default=None,
)
@click.option(
    "--where",
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
    "--dry-run/--no-dry-run",
    help="Perform a trial run with no scan taking place.",
    default=False,
)
def st_sscan(
    ctx,
    trait,
    genotype,
    env,
    covariates,
    kinship,
    where,
    filter_missing,
    filter_maf,
    impute,
    normalize,
    output_dir,
    verbose,
    dry_run,
):
    from os import makedirs
    from os.path import abspath, exists
    from limix._display import session_block, banner, session_line, indent
    from limix.qtl import st_sscan
    from limix.io import fetch
    from .pipeline import Pipeline
    from limix._data import conform_dataset
    from .preprocess import impute as impute_func
    from .preprocess import normalize as normalize_func
    from .preprocess import where as where_func

    print(banner())

    if ctx.obj is None:
        ctx.obj = {"preprocess": []}

    output_dir = abspath(output_dir)
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
            i for i in oparams if i[0] in ["impute", "normalize", "where"]
        ]

        for p in preproc_params:
            if p[0] == "where":
                pipeline.append(where_func, p[1])
            elif p[0] == "normalize":
                pipeline.append(normalize_func, p[1])
            elif p[0] == "impute":
                pipeline.append(impute_func, p[1])

        pipeline.run()

    if dry_run:
        print("Exiting early because of dry-run option.")
        return

    # if "K" not in data:
    #     data["K"] = None
    # try:
    st_sscan(data["G"], data["y"], E=data["E"], M=None, tests=None, verbose=verbose)
    # except Exception as e:
    #     print_exc(traceback.format_stack(), e)
    #     sys.exit(1)

    # with session_line("Saving results to `{}`... ".format(output_dir)):
    #     model.to_csv(join(output_dir, "null.csv"), join(output_dir, "alt.csv"))
    # pass


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
