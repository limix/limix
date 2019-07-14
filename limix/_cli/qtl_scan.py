import click
from click import Path
from loguru import logger

from ._click import limix_command
from ._input import QTLInputData
from ._misc import verbose_option


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
@click.option(
    "--trait", help="Trait file.", default=None, type=Path(exists=True, dir_okay=False)
)
@click.option(
    "--trait-name",
    help="Comma-separated trait names to be used. It defaults to use all traits.",
    default=None,
)
@click.option("--bfile", help="BED/FAM/BIM files prefix.", default=None)
@click.option(
    "--bed", help="BED file.", default=None, type=Path(exists=True, dir_okay=False)
)
@click.option(
    "--fam", help="FAM file.", default=None, type=Path(exists=True, dir_okay=False)
)
@click.option(
    "--bim", help="BIM file.", default=None, type=Path(exists=True, dir_okay=False)
)
@click.option(
    "--grm", help="GRM file.", default=None, type=Path(exists=True, dir_okay=False)
)
@click.option(
    "--rel", help="REL file.", default=None, type=Path(exists=True, dir_okay=False)
)
@click.option(
    "--outdir",
    help="Specify the output directory path.",
    default="output",
    type=Path(file_okay=False, writable=True),
)
@verbose_option
def scan(
    ctx, method, trait, trait_name, bfile, bed, fam, bim, grm, rel, outdir, verbose
):
    import pathlib
    from limix._display import session_line, session_block
    from ._toml import write_limix_toml
    from ._misc import context_info, setup_outdir, setup_logger

    outdir = setup_outdir(outdir)
    setup_logger(outdir)
    write_limix_toml(outdir)

    print(context_info())

    with session_block("Input reading"):
        input = QTLInputData()
        input.set_opt("trait", filepath=trait)
        input.set_opt("trait-name", trait_name=trait_name)
        input.set_opt("bfile", bfile_prefix=bfile)
        input.set_opt("bed", bed_filepath=bed, fam_filepath=fam, bim_filepath=bim)
        input.set_opt("grm", filepath=grm)
        input.set_opt("rel", filepath=rel)
        input.set_opt("outdir", outdir=outdir)

        if input.traits is None:
            raise click.UsageError("no trait has been specified.")

        if input.genotype is None:
            raise click.UsageError("no variant has been specified.")

    print(input)

    for i, task in enumerate(input.tasks(method)):
        task.run(input, i, verbose)
