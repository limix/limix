from optparse import OptionParser


def entry_point():
    from ..mtset.core.postprocessCore import postprocess

    parser = OptionParser()
    parser.add_option("--resdir", dest='resdir', type=str, default='./')
    parser.add_option("--outfile", dest='outfile', type=str, default=None)
    parser.add_option("--tol", dest='tol', type=float, default=4e-3)
    (options, args) = parser.parse_args()
    postprocess(options)
