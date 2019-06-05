from functools import reduce

import click


def limix_command(jointly=[], either=[]):
    class CommandOptionsTogether(click.Command):
        def invoke(self, ctx):
            joint = [list(t) for t in jointly]
            eit = [list(t) for t in either]

            for opts in joint:
                if not _all_true_or_false([ctx.params[opt] is None for opt in opts]):
                    opts = ", ".join([f"--{o}" for o in opts])
                    msg = f"The options [{opts}] must be jointly provided."
                    raise click.ClickException(msg)

            for opts in eit:
                if sum([ctx.params[opt] is not None for opt in opts]) > 1:
                    opts = ", ".join([f"--{o}" for o in opts])
                    msg = f"The options [{opts}] are mutually exclusive."
                    raise click.ClickException(msg)

            super(CommandOptionsTogether, self).invoke(ctx)

    return CommandOptionsTogether


def _all_true_or_false(vals):
    vals = [bool(v) for v in vals]
    return not (reduce(lambda a, b: a & b, vals) ^ reduce(lambda a, b: a | b, vals))
