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
    import os
    from limix._display import session_text, banner, timer_text

    print(banner())
    output_dir = os.path.abspath(output_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    pheno_fetch = limix.io.get_fetch_spec(phenotypes_file)
    geno_fetch = limix.io.get_fetch_spec(genotype_file)

    if verbose:
        print("Phenotype file type: {}".format(pheno_fetch["filetype"]))
        print("Genotype file type: {}".format(geno_fetch["filetype"]))

    y = limix.io.fetch("trait", pheno_fetch, verbose=verbose)
    if verbose:
        print("\n{}\n".format(y))

    G = limix.io.fetch("genotype", geno_fetch, verbose=verbose)
    if verbose:
        print("\n{}\n".format(G))

    data = {"y": y, "G": G}

    with session_text("preprocessing", disable=not verbose):
        _preprocessing(data, filter, filter_missing, filter_maf, impute, verbose)

    try:
        model = limix.qtl.scan(data["G"], data["y"], lik, verbose=verbose)
    except Exception as e:
        from limix import _exception

        _exception.print_exc(traceback.format_stack(), e)
        sys.exit(1)

    with timer_text("Saving results to `{}`... ".format(output_dir)):
        model.to_csv(
            os.path.join(output_dir, "null.csv"), os.path.join(output_dir, "alt.csv")
        )


def _preprocessing(data, filter, filter_missing, filter_maf, impute, verbose):
    from limix._data import conform_dataset
    from .._display import timer_text

    layout = _LayoutChange()

    for target in data.keys():
        layout.append(target, "initial", data[target].shape)

    with timer_text("Matching samples... "):
        data = conform_dataset(data["y"], G=data["G"])
    data = {k: v for k, v in data.items() if v is not None}

    for target in data.keys():
        layout.append(target, "sample match", data[target].shape)

    if data["y"].sample.size == 0:
        print(layout.to_string())
        raise RuntimeError(
            "Exiting early because there is no sample left."
            + " Please, check your sample ids."
        )

    for i, f in enumerate(filter):
        _process_filter(f, data)
        for target in data.keys():
            layout.append(target, "filter {}".format(i), data[target].shape)
            if data["y"].sample.size == 0:
                print(layout.to_string())
                raise RuntimeError("Exiting early because there is no sample left.")

    for f in filter_missing:
        with timer_text("Applying `{}`... ".format(f)):
            _process_filter_missing(f, data)
            if data["y"].sample.size == 0:
                print(layout.to_string())
                raise RuntimeError("Exiting early because there is no sample left.")

    if filter_maf is not None:
        with timer_text("Removing candidates with MAF<{}... ".format(filter_maf)):
            data["G"] = _process_filter_maf(float(filter_maf), data["G"])

        for target in data.keys():
            layout.append(target, "maf filter", data[target].shape)

        if data["G"].candidate.size == 0:
            print(layout.to_string())
            raise RuntimeError("Exiting early because there is no candidate left.")

    for imp in impute:
        _process_impute(imp, data)

    print(layout.to_string())


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


def _process_filter_maf(maf, G):
    import limix

    mafs = limix.qc.compute_maf(G)
    ok = mafs >= maf
    return G.isel(candidate=ok)


def _process_impute(expr, data):
    breakpoint()
    elems = [e.strip() for e in expr.strip().split(":")]
    if len(elems) < 2 or len(elems) > 3:
        raise ValueError("Missing filter syntax error.")

    target = short_name(elems[0])
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


def _print_data_stats(data):
    for k in sorted(data.keys()):
        print(k)
        print(data[k].shape)


class _LayoutChange(object):
    def __init__(self):
        self._targets = {}
        self._steps = ["sentinel"]

    def append(self, target, step, shape):
        if target not in self._targets:
            self._targets[target] = {}

        self._targets[target][step] = shape
        if step != self._steps[-1]:
            self._steps.append(step)

    def to_string(self):
        from texttable import Texttable

        table = Texttable()
        header = [""]
        shapes = {k: [k] for k in self._targets.keys()}

        for step in self._steps[1:]:
            header.append(step)
            for target in self._targets.keys():
                v = str(self._targets[target].get(step, "n/a"))
                shapes[target].append(v)

        table.header(header)

        table.set_cols_dtype(["t"] * len(header))
        table.set_cols_align(["l"] * len(header))
        table.set_deco(Texttable.HEADER)

        for target in self._targets.keys():
            table.add_row(shapes[target])

        msg = table.draw()

        msg = self._add_caption(msg, "-", "Table: Data layout transformation.")
        return msg + "\n"

    def _add_caption(self, msg, c, caption):
        n = len(msg.split("\n")[-1])
        msg += "\n" + (c * n)
        msg += "\n" + caption
        return msg
