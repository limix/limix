from argparse import ArgumentParser

from limix.io import possible_file_types


def _do_see(filepath, filetype, options):
    import limix

    if filetype == 'hdf5':
        limix.io.hdf5.see_hdf5(
            filepath, show_chunks=options.show_chunks, verbose=options.verbose)
    elif filetype == 'csv':
        limix.io.csv.see(
            filepath, verbose=options.verbose, header=options.header == 'yes')
    elif filetype == 'grm.raw':
        limix.io.plink.see_kinship(filepath, verbose=options.verbose)
    elif filetype == 'bed':
        limix.io.plink.see_bed(filepath, verbose=options.verbose)
    elif filetype == 'image':
        limix.plot.see_image(filepath, verbose=options.verbose)
    else:
        print("Unknown file type: %s" % filepath)


def do_see(args):
    import limix

    if args.type is None:
        ft = limix.io.file_type(args.file)
    else:
        ft = args.type

    try:
        _do_see(args.file, ft, args)
    except FileNotFoundError as e:
        print(e)


def see_parser(parser):
    parser.add_argument('file', help='file path')
    parser.add_argument(
        '--quiet', '-q', help='quiet', dest='verbose', action='store_false')
    parser.add_argument(
        '--verbose', help='verbose', dest='verbose', action='store_true')
    parser.add_argument(
        '--show-chunks',
        dest='show_chunks',
        action='store_true',
        help='show chunk information for hdf5 files')
    parser.add_argument(
        '--header',
        choices=['yes', 'no'],
        default='yes',
        help='parse header from file')

    msg = 'specify file type: %s' % ', '.join(possible_file_types())
    parser.add_argument('--type', dest='type', help=msg)

    parser.set_defaults(show_chunks=False, header='yes', verbose=False)
    parser.set_defaults(func=do_see)
    return parser


def do_download(args):
    import limix

    limix.util.download(args.url, verbose=args.verbose)


def download_parser(parser):
    parser.add_argument('url', help='url path')
    parser.add_argument('--quiet', '-q', help='quiet', action='store_true')

    parser.set_defaults(verbose=True)
    parser.set_defaults(func=do_download)
    return parser


def do_extract(args):
    import limix

    limix.util.extract(args.filepath, verbose=args.verbose)


def extract_parser(parser):
    parser.add_argument('filepath', help='file path')
    parser.add_argument('--quiet', '-q', help='quiet', action='store_true')

    parser.set_defaults(verbose=True)
    parser.set_defaults(func=do_extract)
    return parser


def main(args=None):
    try:
        entry(args)
    except SystemExit:
        pass


def entry(args=None):
    p = ArgumentParser(prog="limix")

    subparsers = p.add_subparsers(title='subcommands')
    see_parser(subparsers.add_parser('see'))
    download_parser(subparsers.add_parser('download'))
    extract_parser(subparsers.add_parser('extract'))

    args = p.parse_args(args=args)
    if hasattr(args, 'func'):
        func = args.func
        del args.func
        func(args)
    else:
        p.print_help()
