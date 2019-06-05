import click
from click import Path

from ._click import limix_command
from ._input import InputData
from ._misc import OrderedCommand, ordered_params, verbose_option


@click.command(cls=limix_command([("bed", "fam", "bim")]))
@click.pass_context
@click.option("--pheno", help="Phenotype file.", default=None, type=Path(exists=True))
@click.option("--bfile", help="BED/FAM/BIM files prefix.", default=None)
@click.option("--bed", help="BED file.", default=None, type=Path(exists=True))
@click.option("--fam", help="FAM file.", default=None, type=Path(exists=True))
@click.option("--bim", help="BIM file.", default=None, type=Path(exists=True))
@click.option("--grm", help="GRM file.", default=None, type=Path(exists=True))
@click.option("--rel", help="REL file.", default=None, type=Path(exists=True))
@verbose_option
# @click.option(
#     "--dry-run/--no-dry-run",
#     help="Perform a trial run with no scan taking place.",
#     default=False,
# )
def scan(ctx, pheno, bfile, bed, fam, bim, grm, rel, verbose):
    # if ctx.obj is None:
    #     ctx.obj = {"preprocess": []}

    # params = ordered_params(ctx)
    p = InputData()
    if pheno is not None:
        p.set_pheno(pheno)

    if bfile is not None:
        p.set_bfile(bfile)

    if bed is not None:
        p.set_bed(bed, fam, bim)

    if grm is not None:
        p.set_grm(grm)

    if rel is not None:
        p.set_rel(rel)

    # print("Covariates")
    # print(p.input_data.covariates)

    # print("Genotype")
    # print(p.input_data.genotype)

    # print("Kinship")
    # print(p.input_data.kinship)

    # print("Trait")
    # print(p.input_data.phenotypes)

    if p.phenotypes is None:
        raise click.UsageError("No phenotype has been specified.")

    if p.genotype is None:
        raise click.UsageError("No variant has been specified.")

    # for y in p.input_data.phenotypes:
    #     print("Trait {}".format(y.name))
    #     data = {
    #         "y": y,
    #         "M": p.input_data.covariates,
    #         "K": p.input_data.kinship,
    #         "G": p.input_data.genotype,
    #     }
    #     with session_line("Matching samples... "):
    #         data = conform_dataset(**data)
    #     data = {k: v for k, v in data.items() if v is not None}
