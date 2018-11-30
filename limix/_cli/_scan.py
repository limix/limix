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
    from os import makedirs
    from os.path import abspath, exists, join
    import traceback
    from .._display import session_block, banner, session_line, print_exc
    from ._spec import parse_fetch_spec

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

    with session_block("preprocessing", disable=not verbose):
        data = _preprocessing(data, filter, filter_missing, filter_maf, impute, verbose)

    try:
        model = limix.qtl.st_scan(
            data["G"], data["y"], lik, K=data["K"], verbose=verbose
        )
    except Exception as e:
        print_exc(traceback.format_stack(), e)
        sys.exit(1)

    with session_line("Saving results to `{}`... ".format(output_dir)):
        model.to_csv(join(output_dir, "null.csv"), join(output_dir, "alt.csv"))


def _preprocessing(data, filter, filter_missing, filter_maf, impute, verbose):
    from limix._data import conform_dataset
    from .._display import session_line

    layout = _LayoutChange()

    for target in data.keys():
        layout.append(target, "initial", data[target].shape)

    with session_line("Matching samples... "):
        data = conform_dataset(**data)
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
        data = _process_filter(f, data)
        for target in data.keys():
            layout.append(target, "filter {}".format(i), data[target].shape)
            if data["y"].sample.size == 0:
                print(layout.to_string())
                raise RuntimeError("Exiting early because there is no sample left.")

    for f in filter_missing:
        with session_line("Applying `{}`... ".format(f)):
            _process_filter_missing(f, data)
            if data["y"].sample.size == 0:
                print(layout.to_string())
                raise RuntimeError("Exiting early because there is no sample left.")

    if filter_maf is not None:
        with session_line("Removing candidates with MAF<{}... ".format(filter_maf)):
            data["G"] = _process_filter_maf(float(filter_maf), data["G"])

        for target in data.keys():
            layout.append(target, "maf filter", data[target].shape)

        if data["G"].candidate.size == 0:
            print(layout.to_string())
            raise RuntimeError("Exiting early because there is no candidate left.")

    for imp in impute:
        with session_line("Imputting missing values (`{}`)... ".format(imp)):
            data = _process_impute(imp, data)

    print(layout.to_string())

    return data


def _process_filter(expr, data):
    from .._bits.xarray import query
    from .._data import to_short_data_name

    elems = [e.strip() for e in expr.strip().split(":")]
    if len(elems) < 2 or len(elems) > 3:
        raise ValueError("Filter syntax error.")
    target = elems[0]
    expr = elems[1]
    n = to_short_data_name(target)
    data[n] = query(data[n], expr)
    return data


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
    from limix import compute_maf

    mafs = compute_maf(G)
    ok = mafs >= maf
    return G.isel(candidate=ok)


def _process_impute(expr, data):
    from .._data import to_short_data_name, dim_hint_to_name, dim_name_to_hint

    elems = [e.strip() for e in expr.strip().split(":")]
    if len(elems) < 2 or len(elems) > 3:
        raise ValueError("Missing filter syntax error.")

    target = to_short_data_name(elems[0])
    dim = elems[1]

    if len(elems) == 3:
        method = elems[2]
    else:
        method = "mean"

    def in_dim(X, dim):
        return dim_hint_to_name(dim) in X.dims or dim_name_to_hint(dim) in X.dims

    X = data[target]
    if not in_dim(X, dim):
        raise ValueError("Unrecognized dimension: {}.".format(dim))

    if method == "mean":
        axis = next(i for i in range(len(X.dims)) if in_dim(X, dim))
        if axis == 0:
            X = limix.qc.impute.mean_impute(X.T).T
        else:
            X = limix.qc.impute.mean_impute(X)
    else:
        raise ValueError("Unrecognized imputation method: {}.".format(method))

    data[target] = X

    return data


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
        return msg

    def _add_caption(self, msg, c, caption):
        n = len(msg.split("\n")[-1])
        msg += "\n" + (c * n)
        msg += "\n" + caption
        return msg
