import traceback
import sys
import limix
import click
import os


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
    multiple=True,
)
@click.option(
    "--impute",
    help=("Impute missing values for phenotype, genotype, and covariate."),
    multiple=True,
)
@click.option(
    "--output-dir", help="Specify the output directory path.", default="output"
)
@click.option(
    "--verbose/--quiet", "-v/-q", help="Enable or disable verbose mode.", default=True
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
    output_dir,
    verbose,
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

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    pheno_fetch = limix.io.get_fetch_specification(phenotypes_file)
    geno_fetch = limix.io.get_fetch_specification(genotype_file)

    y = limix.io.fetch_phenotype(pheno_fetch, verbose=verbose)
    G = limix.io.fetch_genotype(geno_fetch, verbose=verbose)

    data = {"phenotype": y, "genotype": G}

    for f in filter:
        _process_filter(f, data)

    for f in filter_missing:
        _process_filter_missing(f, data)

    for f in filter_maf:
        _process_filter_maf(f, data)

    for imp in impute:
        _process_impute(imp, data)

    try:
        model = limix.qtl.scan(
            data["genotype"], data["phenotype"], lik, verbose=verbose
        )
    except Exception as e:
        limix._exception.print_exc(traceback.format_stack(), e)
        sys.exit(1)

    if verbose:
        print(model)

    model.to_csv(
        os.path.join(output_dir, "null.csv"), os.path.join(output_dir, "alt.csv")
    )


def _process_filter(expr, data):
    elems = [e.strip() for e in expr.strip().split(":")]
    if len(elems) < 2 or len(elems) > 3:
        raise ValueError("Filter syntax error.")


def _process_filter_missing(expr, data):
    elems = [e.strip() for e in expr.strip().split(":")]
    if len(elems) < 2 or len(elems) > 3:
        raise ValueError("Missing filter syntax error.")

    target = elems[0]
    dim = elems[1]

    if len(elems) == 3:
        how = elems[2]
    else:
        how = "any"

    data[target] = data[target].dropna(dim, how)


def _process_filter_maf(expr, data):
    maf = float(expr)
    G = data["genotype"]
    data["genotype"] = G


def _process_impute(expr, data):
    elems = [e.strip() for e in expr.strip().split(":")]
    if len(elems) < 2 or len(elems) > 3:
        raise ValueError("Missing filter syntax error.")

    target = elems[0]
    dim = elems[1]

    if len(elems) == 3:
        method = elems[2]
    else:
        method = "mean"

    X = data[target]
    if dim not in X.dims:
        raise ValueError("Unrecognized dimension: {}.".format(dim))

    if method == "mean":
        if X.dims[0] == dim:
            X = limix.qc.impute.mean_impute(X.T).T
        else:
            X = limix.qc.impute.mean_impute(X)
    else:
        raise ValueError("Unrecognized imputation method: {}.".format(method))

    data[target] = X
