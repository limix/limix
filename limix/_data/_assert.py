def assert_target(target):
    from ._conf import CONF

    if target not in CONF["targets"]:
        raise ValueError(f"Uknown target `{target}`.")


def assert_filetype(filetype):
    from ._conf import CONF

    if filetype not in CONF["filetypes"]:
        raise ValueError(f"Uknown filetype `{filetype}`.")
