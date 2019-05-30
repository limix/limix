import click

from ._misc import OrderedCommand, verbose_option, ordered_params
from ._input import ProcessInputData


@click.command(cls=OrderedCommand)
@click.pass_context
@click.option("--bed", help="BED file.", default=None)
@click.option("--fam", help="FAM file.", default=None)
@click.option("--bim", help="BIM file.", default=None)
@click.option("--grm", help="GRM file.", default=None)
@click.option("--rel", help="REL file.", default=None)
@click.option("--pheno", help="Phenotype file.", default=None)
@click.option("--bfile", help="BED/FAM/BIM files prefix.", default=None)
@verbose_option
# @click.option(
#     "--dry-run/--no-dry-run",
#     help="Perform a trial run with no scan taking place.",
#     default=False,
# )
def scan(ctx, **_):
    from limix._display import session_line
    from limix._data import conform_dataset

    if ctx.obj is None:
        ctx.obj = {"preprocess": []}

    params = ordered_params(ctx)
    p = ProcessInputData(params)

    print("Covariates")
    print(p.input_data.covariates)

    print("Genotype")
    print(p.input_data.genotype)

    print("Kinship")
    print(p.input_data.kinship)

    print("Trait")
    print(p.input_data.phenotypes)

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
