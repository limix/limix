from optparse import OptionParser

from ..mtset.core.definesetsCore import definesets


def entry_point():
    parser = OptionParser()
    parser.add_option("--bfile", dest='bfile', type=str, default=None)
    parser.add_option("--wfile", dest='wfile', type=str, default=None)

    parser.add_option(
        "--sliding_window",
        action="store_true",
        dest="sliding_wind",
        default=False)

    parser.add_option(
        "--filter_sets",
        action="store_true",
        dest='filter_sets',
        default=False)

    parser.add_option(
        "--window_size", dest='window_size', type=int, default=3e4)

    parser.add_option("--step", dest='step', type=int, default=15e3)

    parser.add_option("--minSnps", dest='minSnps', type=int, default=1)
    parser.add_option("--maxSnps", dest='maxSnps', type=int, default=None)
    parser.add_option("--chrom", dest='chrom', type=int, default=None)

    parser.add_option("--iwfile", dest='iwfile', type=str, default=None)

    (options, args) = parser.parse_args()

    definesets(options)
