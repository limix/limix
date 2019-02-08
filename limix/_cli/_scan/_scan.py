import sys

import click
import limix


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
        " to use in the analysis."
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
        "Drop out candidates having a minor allele frequency below the provided threshold."
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
    from os.path import abspath, exists, join
    import traceback
    from limix._display import session_block, banner, session_line, print_exc
    from .._spec import parse_fetch_spec
    from .._preprocess import _preprocessing

    print(banner())

    output_dir = abspath(output_dir)
    if not exists(output_dir):
        makedirs(output_dir)

    fetch = {
        "phenotype": parse_fetch_spec(phenotypes_file),
        "genotype": parse_fetch_spec(genotype_file),
    }
    if kinship_file is not None:
        fetch["kinship"] = parse_fetch_spec(kinship_file)

    if verbose:
        print("Phenotype file type: {}".format(fetch["phenotype"]["filetype"]))
        print("Genotype file type: {}".format(fetch["genotype"]["filetype"]))
        if "kinship" in fetch:
            print("Kinship file type: {}".format(fetch["kinship"]["filetype"]))

    y = limix.io.fetch("trait", fetch["phenotype"], verbose=verbose)
    y = _fix_trait_dims(y)
    if verbose:
        print("\n{}\n".format(y))

    G = limix.io.fetch("genotype", fetch["genotype"], verbose=verbose)
    if verbose:
        print("\n{}\n".format(G))

    data = {"y": y, "G": G, "K": None}
    if kinship_file is not None:
        K = limix.io.fetch("covariance", fetch["kinship"], verbose=verbose)
        if verbose:
            print("\n{}\n".format(K))
        data["K"] = K
    if data["K"] is None:
        del data["K"]

    with session_block("preprocessing", disable=not verbose):
        data = _preprocessing(
            data, filter, filter_missing, filter_maf, impute, normalize, verbose
        )

    if "K" not in data:
        data["K"] = None
    try:
        model = limix.qtl.st_scan(
            data["G"], data["y"], lik, K=data["K"], verbose=verbose
        )
    except Exception as e:
        print_exc(traceback.format_stack(), e)
        sys.exit(1)

    with session_line("Saving results to `{}`... ".format(output_dir)):
        model.to_csv(join(output_dir, "null.csv"), join(output_dir, "alt.csv"))


def _fix_trait_dims(y):
    if y.ndim == 1:
        if y.name != "trait":
            raise RuntimeError(
                "The name of an unidimensional trait array should be" " 'trait'."
            )
        if y.dims[0] != "sample":
            y = y.rename({y.dims[0]: "sample"})
    return y
