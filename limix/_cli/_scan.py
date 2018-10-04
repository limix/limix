import traceback
import sys
import limix
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
    from ._exception import print_exc

    pheno_filepath, pheno_type = limix.io.detect_file_type(phenotypes_file)

    geno_filepath, geno_type = limix.io.detect_file_type(genotype_file)

    y = limix.io.fetch_phenotype(pheno_filepath, pheno_type)
    G = limix.io.fetch_genotype(geno_filepath, geno_type)

    data = {"phenotype": y, "genotype": G}

    for filt in filter:
        target = _get_filter_target(filt)
        data[target] = _dispath_process_filter[target](data[target], filt)

    try:
        r = limix.qtl.scan(data["genotype"], data["phenotype"], lik)
    except Exception as e:
        print_exc(traceback.format_stack(), e)
        sys.exit(1)
    print(r)


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
