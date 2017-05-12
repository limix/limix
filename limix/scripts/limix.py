from argparse import ArgumentParser


def do_see(args):
    import limix

    ft = limix.io.file_type(args.file)

    if ft == 'hdf5':
        limix.io.hdf5.see(args.file, show_chunks=args.show_chunks)
    else:
        print("Unknown file type: %s." % args.file)


def parse_see(args):
    p = ArgumentParser(prog='limix')
    p.add_argument('file')
    p.add_argument(['-h', '--help'], dest='help', action='store_true')
    p.add_argument('--show-chunks', dest='show_chunks', action='store_true')
    p.set_defaults(help=False)
    p.set_defaults(show_chunks=False)

    args = p.parse_args(args)

    if args.help:
        p.print_help()
    else:
        do_see(args)


def entry_point():
    p = ArgumentParser()

    sub = p.add_subparsers(title='subcommands')

    s = sub.add_parser('see')
    s.add_argument('file', help='file path')
    s.add_argument('--show-chunks', dest='show_chunks', action='store_true',
                   help='show chunk information')
    s.set_defaults(show_chunks=False)

    # s.set_defaults(func=parse_see)
    #
    args = p.parse_args()
    # args, rargs = p.parse_known_args()

    #
    # if hasattr(args, 'func'):
    #     func = args.func
    #     del args.func
    #     func(rargs)
    # else:
    #     p.print_help()
