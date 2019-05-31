from functools import reduce

import click


def limix_command(jointly=[]):
    class CommandOptionsTogether(click.Command):
        def invoke(self, ctx):
            together = [list(t) for t in jointly]

            for opts in together:
                if not _all_true_or_false([ctx.params[opt] is None for opt in opts]):
                    opts = ", ".join([f"--{o}" for o in opts])
                    msg = f"The options [{opts}] must be jointly provided."
                    raise click.ClickException(msg)

            super(CommandOptionsTogether, self).invoke(ctx)

    return CommandOptionsTogether


def _all_true_or_false(vals):
    vals = [bool(v) for v in vals]
    return not (reduce(lambda a, b: a & b, vals) ^ reduce(lambda a, b: a | b, vals))
