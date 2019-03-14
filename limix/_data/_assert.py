def assert_target(target):
    from ._conf import CONF

    if target not in CONF["targets"]:
        raise ValueError(f"Unknown target `{target}`.")


def assert_filetype(filetype):
    from ._conf import CONF

    if filetype not in CONF["filetypes"]:
        raise ValueError(f"Unknown filetype `{filetype}`.")


def assert_likelihood(likname):
    from ._conf import CONF

    if likname not in CONF["likelihoods"]:
        msg = "Unrecognized likelihood name: {}.\n".format(likname)
        msg += "Valid names are: {}.".format(list(CONF["likelihoods"]))
        raise ValueError(msg)
