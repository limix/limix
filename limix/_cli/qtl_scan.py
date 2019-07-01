import click
from click import Path

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
# @click.option(
#     "--dry-run/--no-dry-run",
#     help="Perform a trial run with no scan taking place.",
#     default=False,
# )
def scan(
    ctx, method, trait, trait_name, bfile, bed, fam, bim, grm, rel, outdir, verbose
):
    import time
    import os
    import sys
    from pathlib import Path
    from os.path import join, abspath, exists
    from os import makedirs
    from limix._display import banner, session_line, session_block

    def _curdate():
        return time.strftime('%l:%M:%S%p %Z on %b %d, %Y')

    print(banner())
    start_date = _curdate()
    outdir = Path(outdir)
    workdir = os.getcwd()

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

    _single_trait(input, verbose)

    end_date = _curdate()
    cmdline = " ".join(sys.argv)
    info = f"workdir={workdir}\ncmdline={cmdline}\n"
    info += f"start_date={start_date}\nend_date={end_date}"
    with open(input.outdir / "info.txt", "w") as f:
        f.write(info)


def _single_trait(input, verbose):
    import limix
    from limix._data import conform_dataset
    from limix._display import session_line

    Y = input.traits.T
    covariates = input.covariates
    kinship = input.kinship
    genotype = input.genotype

    for y in Y:
        data = {"y": y, "M": covariates, "K": kinship, "G": genotype}
        data = conform_dataset(**data)
        data = {k: v for k, v in data.items() if v is not None}

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
        outdir = input.outdir / y.trait.values.item()
        if not outdir.exists():
            outdir.mkdir()

        outdir_stem = str(outdir.absolute())
        with session_line(f"Saving results to `{outdir_stem}`... "):
            res.to_csv(
                outdir / "h0_effsizes.csv",
                outdir / "h0_variances.csv",
                outdir / "h2_effsizes.csv",
                outdir / "stats.csv",
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
