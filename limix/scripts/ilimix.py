import cmd
from . import limix as limix_cmd
import shlex

class ILimix(cmd.Cmd):

    def do_see(self, cmdline):
        limix_cmd.parse_see(shlex.split(cmdline))

    def do_EOF(self, _):
        return True

    def do_exit(self, *_):
        return True

def entry_point():
    ILimix().cmdloop()

if __name__ == '__main__':
    entry_point()
