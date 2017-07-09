from optparse import OptionParser


def entry_point():
    from ..mtSet.core.preprocessCore import preprocess

    parser = OptionParser()
    parser.add_option("--bfile", dest='bfile', type=str, default=None)
    parser.add_option("--cfile", dest='cfile', type=str, default=None)
    parser.add_option("--pfile", dest='pfile', type=str, default=None)
    parser.add_option("--nfile", dest='nfile', type=str, default=None)
    parser.add_option("--ffile", dest='ffile', type=str, default=None)
    parser.add_option("--trait_idx", dest='trait_idx', type=str, default=None)

    parser.add_option(
        "--compute_covariance",
        action="store_true",
        dest="compute_cov",
        default=False)
    parser.add_option("--compute_PCs", dest="compute_PCs", default=0, type=int)

    parser.add_option(
        "--plink_path", dest='plink_path', type=str, default='plink')
    parser.add_option("--sim_type", dest='sim_type', type=str, default='RRM')

    parser.add_option(
        "--fit_null", action="store_true", dest="fit_null", default=False)

    (options, args) = parser.parse_args()

    preprocess(options)
