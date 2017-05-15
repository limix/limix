import cmd
import shlex
from argparse import ArgumentParser

from . import limix as limix_cmd


class ILimix(cmd.Cmd):
    def do_see(self, cmdline):
        p = limix_cmd.see_parser(ArgumentParser())
        args = p.parse_args(shlex.split(cmdline))

        limix_cmd.do_see(args)

    def do_EOF(self, _):
        return True

    def do_exit(self, *_):
        return True


def entry_point():
    ILimix().cmdloop()


if __name__ == '__main__':
    entry_point()
