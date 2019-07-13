import sys


def running_environment():
    import os
    from time import strftime
    from ._aligned import AlignedText
    from limix import __version__

    pyver = sys.version.split("\n")[0].strip()
    workdir = os.getcwd()
    start_date = strftime('%I:%M:%S%p %Z on %b %d, %Y')
    cmdline = " ".join(sys.argv)

    aligned = AlignedText(" ")
    aligned.add_item("Limix", __version__)
    aligned.add_item("Python", pyver)
    aligned.add_item("Date", start_date)
    aligned.add_item("Workdir", workdir)
    aligned.add_item("Cmdline", cmdline)
    return aligned.draw()
