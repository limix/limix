from argparse import ArgumentParser

from limix.io import possible_file_types


def do_see(args):
    import limix

    if args.type is None:
        ft = limix.io.file_type(args.file)
    else:
        ft = args.type

    if ft == 'hdf5':
        limix.io.hdf5.see_hdf5(
            args.file, show_chunks=args.show_chunks, verbose=not args.quiet)
    elif ft == 'csv':
        limix.io.csv.see(args.file, verbose=not args.quiet)
    elif ft == 'grm.raw':
        limix.io.plink.see_kinship(args.file, verbose=not args.quiet)
    elif ft == 'bed':
        limix.io.plink.see_bed(args.file, verbose=not args.quiet)
    elif ft == 'image':
        limix.plot.see_image(args.file, verbose=not args.quiet)
    else:
        print("Unknown file type: %s" % args.file)


def see_parser(parser):
    parser.add_argument('file', help='file path')
    parser.add_argument('--quiet', '-q', help='quiet', action='store_true')
    parser.add_argument(
        '--show-chunks',
        dest='show_chunks',
        action='store_true',
        help='show chunk information for hdf5 files')

    msg = 'specify file type: %s' % ', '.join(possible_file_types())
    parser.add_argument('--type', dest='type', help=msg)
    parser.set_defaults(show_chunks=False, verbose=True)
    parser.set_defaults(func=do_see)
    return parser


def do_download(args):
    import limix

    limix.util.download(args.url, verbose=not args.quiet)


def download_parser(parser):
    parser.add_argument('url', help='url path')
    parser.add_argument('--quiet', '-q', help='quiet', action='store_true')

    parser.set_defaults(verbose=True)
    parser.set_defaults(func=do_download)
    return parser


def do_extract(args):
    import limix

    limix.util.extract(args.filepath, verbose=not args.quiet)


def extract_parser(parser):
    parser.add_argument('filepath', help='file path')
    parser.add_argument('--quiet', '-q', help='quiet', action='store_true')

    parser.set_defaults(verbose=True)
    parser.set_defaults(func=do_extract)
    return parser


def entry_point():
    p = ArgumentParser()

    subparsers = p.add_subparsers(title='subcommands')
    see_parser(subparsers.add_parser('see'))
    download_parser(subparsers.add_parser('download'))
    extract_parser(subparsers.add_parser('extract'))

    args = p.parse_args()
    if hasattr(args, 'func'):
        func = args.func
        del args.func
        func(args)
    else:
        p.print_help()
