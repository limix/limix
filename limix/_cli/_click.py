import click


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
