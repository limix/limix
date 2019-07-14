import click
from loguru import logger


def verbose_option(func):
    return click.option(
        "--verbose/--quiet",
        "-v/-q",
        help="Enable or disable verbose mode.",
        default=True,
    )(func)


def get_version():
    import pkg_resources
    import re
    from os.path import realpath, dirname, join

    if __name__ == "__main__":
        filepath = join(dirname(realpath(__file__)), "..", "__init__.py")
        with open(filepath, "r", encoding="utf8") as f:
            content = f.read()
    else:
        content = pkg_resources.resource_string(__name__.split(".")[0], "__init__.py")
        content = content.decode()

    c = re.compile(r"__version__ *= *('[^']+'|\"[^\"]+\")")
    m = c.search(content)
    if m is None:
        return "unknown"
    return m.groups()[0][1:-1]


class OrderedCommand(click.Command):
    def parse_args(self, ctx, args):
        parser = self.make_parser(ctx)
        opts, args, param_order = parser.parse_args(args=args)

        for param in click.core.iter_params_for_processing(
            param_order, self.get_params(ctx)
        ):
            value, args = param.handle_parse_result(ctx, opts, args)

        if args and not ctx.allow_extra_args and not ctx.resilient_parsing:
            ctx.fail(
                "Got unexpected extra argument%s (%s)"
                % (
                    len(args) != 1 and "s" or "",
                    " ".join(map(click.utils.make_str, args)),
                )
            )

        ctx.args = args
        ctx.param_order = param_order
        return args


def ordered_params(ctx):
    args_seq = []
    for p in ctx.param_order:
        for opt, val in ctx.params.items():
            if p.name == opt:
                if isinstance(val, tuple):
                    v = val[0]
                    val = val[1:]
                else:
                    v = val
                args_seq.append((opt, v))
                if isinstance(val, tuple):
                    if len(val) == 0:
                        del ctx.params[opt]
                    else:
                        ctx.params[opt] = v
                else:
                    del ctx.params[opt]
                break
            pass
    return args_seq


def context_info():
    import sys
    import os
    from time import strftime
    from limix import __version__
    from limix._display import AlignedText

    pyver = sys.version.split("\n")[0].strip()
    workdir = os.getcwd()
    start_date = strftime("%I:%M:%S%p %Z on %b %d, %Y")
    cmdline = " ".join(sys.argv)

    aligned = AlignedText(" ")
    aligned.add_item("Limix", __version__)
    aligned.add_item("Python", pyver)
    aligned.add_item("Date", start_date)
    aligned.add_item("Workdir", workdir)
    aligned.add_item("Cmdline", cmdline)

    msg = aligned.draw()
    logger.info("\n" + msg)
    return msg


def setup_outdir(outdir):
    from pathlib import Path

    outdir = Path(outdir)
    if not outdir.exists():
        outdir.mkdir()
    return outdir


def setup_logger(outdir):
    from loguru import logger

    logger.remove()
    logger.add(outdir / "limix.log", level="INFO")
