import click

from ._download import download
from ._estimate_kinship import estimate_kinship
from ._extract import extract
from ._remove import remove
from ._scan import scan
from ._see import see


def _get_version():
    import pkg_resources
    import re
    from os.path import realpath, dirname, join

    if __name__ == "__main__":
        filepath = join(dirname(realpath(__file__)), "..", "__init__.py")
        with open(filepath, "r") as f:
            content = f.read()
    else:
        content = pkg_resources.resource_string(__name__.split(".")[0], "__init__.py")
        content = content.decode()

    c = re.compile(r"__version__ *= *('[^']+'|\"[^\"]+\")")
    m = c.search(content)
    if m is None:
        return "unknown"
    return m.groups()[0][1:-1]


@click.group(name="limix", context_settings=dict(help_option_names=["-h", "--help"]))
@click.pass_context
@click.version_option(_get_version())
def cli(ctx):
    pass


cli.add_command(see)
cli.add_command(estimate_kinship)
cli.add_command(download)
cli.add_command(extract)
cli.add_command(scan)
cli.add_command(remove)
