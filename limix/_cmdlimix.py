from argparse import ArgumentParser


def _do_see(filepath, filetype, options):
    import limix

    if filetype == "hdf5":
        limix.io.hdf5.see_hdf5(
            filepath, show_chunks=options.show_chunks, verbose=options.verbose
        )
    elif filetype == "csv":
        limix.io.csv.see(
            filepath, verbose=options.verbose, header=options.header == "yes"
        )
    elif filetype == "grm.raw":
        limix.io.plink.see_kinship(filepath, verbose=options.verbose)
    elif filetype == "bed":
        limix.io.plink.see_bed(filepath, verbose=options.verbose)
    elif filetype == "image":
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


def do_estimate_kinship(args):
    import limix

    input_file = getattr(args, "input-file")

    if args.type is None:
        ft = limix.io.file_type(input_file)
    else:
        ft = args.type

    if args.verbose:
        print("Detected file type: {}".format(ft))

    if ft == "bgen":
        G = limix.io.bgen.fetch_dosage(input_file, verbose=args.verbose)
    elif ft == "bed":
        G = limix.io.plink.fetch_dosage(input_file, verbose=args.verbose)
    else:
        print("Unknown file type: %s" % input_file)

    K = limix.stats.linear_kinship(G, verbose=args.verbose)

    if args.output_file is None:
        output_file = input_file + ".npy"
    else:
        output_file = args.output_file

    oft = limix.io.file_type(output_file)

    if oft == "npy":
        limix.io.npy.save_kinship(output_file, K, verbose=args.verbose)
    else:
        print("Unknown output file type: %s" % output_file)


def estimate_kinship_parser(parser):
    from limix.io import possible_file_types

    parser.add_argument("input-file", help="input file path")
    parser.add_argument(
        "--quiet", "-q", help="quiet", dest="verbose", action="store_false"
    )
    parser.add_argument(
        "--verbose", help="verbose", dest="verbose", action="store_true"
    )
    parser.add_argument(
        "--output-file", dest="output_file", default=None, help="output file path"
    )

    msg = "specify input file type: %s" % ", ".join(possible_file_types())
    parser.add_argument("--type", dest="type", help=msg)

    parser.set_defaults(output_file=None, verbose=True)
    parser.set_defaults(func=do_estimate_kinship)
    return parser


def see_parser(parser):
    from limix.io import possible_file_types

    parser.add_argument("file", help="file path")
    parser.add_argument(
        "--quiet", "-q", help="quiet", dest="verbose", action="store_false"
    )
    parser.add_argument(
        "--verbose", help="verbose", dest="verbose", action="store_true"
    )
    parser.add_argument(
        "--show-chunks",
        dest="show_chunks",
        action="store_true",
        help="show chunk information for hdf5 files",
    )
    parser.add_argument(
        "--header", choices=["yes", "no"], default="yes", help="parse header from file"
    )

    msg = "specify file type: %s" % ", ".join(possible_file_types())
    parser.add_argument("--type", dest="type", help=msg)

    parser.set_defaults(show_chunks=False, header="yes", verbose=False)
    parser.set_defaults(func=do_see)
    return parser


def do_download(args):
    import limix

    limix.sh.download(args.url, verbose=args.verbose)


def download_parser(parser):
    parser.add_argument("url", help="url path")
    parser.add_argument("--quiet", "-q", help="quiet", action="store_true")

    parser.set_defaults(verbose=True)
    parser.set_defaults(func=do_download)
    return parser


def remove_parser(parser):
    parser.add_argument("filepath", help="file path")
    parser.add_argument("--quiet", "-q", help="quiet", action="store_true")

    parser.set_defaults(verbose=True)
    parser.set_defaults(func=do_remove)
    return parser


def do_remove(args):
    import limix

    limix.sh.remove(args.filepath)


def do_extract(args):
    import limix

    limix.sh.extract(args.filepath, verbose=args.verbose)


def extract_parser(parser):
    parser.add_argument("filepath", help="file path")
    parser.add_argument("--quiet", "-q", help="quiet", action="store_true")

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

    subparsers = p.add_subparsers(title="subcommands")
    see_parser(subparsers.add_parser("see"))
    download_parser(subparsers.add_parser("download"))
    remove_parser(subparsers.add_parser("remove"))
    extract_parser(subparsers.add_parser("extract"))
    estimate_kinship_parser(subparsers.add_parser("estimate-kinship"))

    args = p.parse_args(args=args)
    if hasattr(args, "func"):
        func = args.func
        del args.func
        func(args)
    else:
        p.print_help()
