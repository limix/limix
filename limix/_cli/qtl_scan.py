import click
from click import Path

from ._click import limix_command
from ._input import InputData
from ._misc import OrderedCommand, ordered_params, verbose_option


@click.command(
    cls=limix_command(jointly=[("bed", "fam", "bim")], either=[("grm", "rel")])
)
@click.pass_context
@click.option("--pheno", help="Phenotype file.", default=None, type=Path(exists=True))
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
def scan(ctx, pheno, bfile, bed, fam, bim, grm, rel, outdir, verbose):
    import limix
    from os.path import join, abspath, exists
    from os import makedirs
    import sys
    from limix._data import conform_dataset
    from limix._display import banner, session_line, print_exc
    import traceback

    print(banner())

    outdir = abspath(outdir)
    if not exists(outdir):
        makedirs(outdir, exist_ok=True)

    p = InputData()
    p.set_opt("pheno", pheno=pheno)
    p.set_opt("bfile", bfile=bfile)
    p.set_opt("bed", bed=bed, fam=fam, bim=bim)
    p.set_opt("grm", grm=grm)
    p.set_opt("rel", rel=rel)

    if p.phenotypes is None:
        raise click.UsageError("No phenotype has been specified.")

    if p.genotype is None:
        raise click.UsageError("No variant has been specified.")

    for y in p.phenotypes.T:
        print(f"Phenotype: {y.name}")
        data = {"y": y, "M": p.covariates, "K": p.kinship, "G": p.genotype}
        data = conform_dataset(**data)
        data = {k: v for k, v in data.items() if v is not None}

        if "K" not in data:
            data["K"] = None
        try:
            res = limix.qtl.scan(
                data["G"],
                data["y"],
                lik="normal",
                K=data["K"],
                M=data["M"],
                verbose=verbose,
            )
        except Exception as e:
            print_exc(traceback.format_stack(), e)
            sys.exit(1)

        with session_line("Saving results to `{}`... ".format(outdir)):
            res.to_csv(
                join(outdir, "h0_effsizes.csv"),
                join(outdir, "h0_variances.csv"),
                join(outdir, "h2_effsizes.csv"),
                join(outdir, "stats.csv"),
            )
