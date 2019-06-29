import click
from click import Path

from ._click import limix_command
from ._input import QTLInputData
from ._misc import OrderedCommand, ordered_params, verbose_option


@click.command(
    cls=limix_command(jointly=[("bed", "fam", "bim")], either=[("grm", "rel")])
)
@click.pass_context
@click.option(
    "--method",
    help="QTL method. Defaults to `st`.",
    default="st",
    type=click.Choice(["st", "mt", "struct"]),
)
@click.option("--pheno", help="Phenotype file.", default=None, type=Path(exists=True))
@click.option(
    "--pheno-name",
    help="Comma-separated phenotype names to be used in the analysis. It defaults to use all loaded phenotypes.",
    default=None,
)
@click.option("--bfile", help="BED/FAM/BIM files prefix.", default=None)
@click.option("--bed", help="BED file.", default=None, type=Path(exists=True))
@click.option("--fam", help="FAM file.", default=None, type=Path(exists=True))
@click.option("--bim", help="BIM file.", default=None, type=Path(exists=True))
@click.option("--grm", help="GRM file.", default=None, type=Path(exists=True))
@click.option("--rel", help="REL file.", default=None, type=Path(exists=True))
@click.option("--outdir", help="Specify the output directory path.", default="output")
@verbose_option
# @click.option(
#     "--dry-run/--no-dry-run",
#     help="Perform a trial run with no scan taking place.",
#     default=False,
# )
def scan(
    ctx, method, pheno, pheno_name, bfile, bed, fam, bim, grm, rel, outdir, verbose
):
    from os.path import join, abspath, exists
    from os import makedirs
    from limix._display import banner, session_line, session_block

    print(banner())

    outdir = abspath(outdir)
    if not exists(outdir):
        makedirs(outdir, exist_ok=True)

    with session_block("Input reading"):
        p = QTLInputData()
        p.set_opt("pheno", filepath=pheno)
        p.set_opt("pheno-name", pheno_name=pheno_name)
        p.set_opt("bfile", bfile_prefix=bfile)
        p.set_opt("bed", bed_filepath=bed, fam_filepath=fam, bim_filepath=bim)
        p.set_opt("grm", filepath=grm)
        p.set_opt("rel", filepath=rel)

        if p.phenotypes is None:
            raise click.UsageError("no phenotype has been specified.")

        if p.genotype is None:
            raise click.UsageError("no variant has been specified.")

    print(p)

    _single_trait(p, verbose)

    # with session_line(f"Saving results to `{outdir}`... "):
    #     res.to_csv(
    #         join(outdir, "h0_effsizes.csv"),
    #         join(outdir, "h0_variances.csv"),
    #         join(outdir, "h2_effsizes.csv"),
    #         join(outdir, "stats.csv"),
    #     )


def _single_trait(input, verbose):
    import limix
    from limix._data import conform_dataset

    Y = input.phenotypes.T
    covariates = input.covariates
    kinship = input.kinship
    genotype = input.genotype

    for y in Y:
        data = {"y": y, "M": covariates, "K": kinship, "G": genotype}
        data = conform_dataset(**data)
        data = {k: v for k, v in data.items() if v is not None}
        # print_trait(y_given, data["y"])
        # return

        if "K" not in data:
            data["K"] = None

        res = limix.qtl.scan(
            data["G"],
            data["y"],
            lik="normal",
            K=data["K"],
            M=data["M"],
            verbose=verbose,
        )


def _multi_trait(input):
    pass


def _struct_lmm(input):
    pass


def print_trait(y_given, y_used):
    from limix._display import summarize_list_repr

    print(f"Phenotype: {y_given.name}")
    n_given = len(y_given)
    samples_given = summarize_list_repr(sorted(y_given.sample.values.tolist()), 5)
    print(f"  Samples given ({n_given}): {samples_given}")
    n_used = len(y_used)
    samples_used = summarize_list_repr(sorted(y_used.sample.values.tolist()), 5)
    print(f"  Samples used ({n_used}): {samples_used}")
    pass
