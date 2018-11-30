from ._conf import CONF


def is_data_name(name):
    return name in CONF["data_names"]


def is_short_data_name(name):
    return name in CONF["short_data_names"]


def to_short_data_name(name):
    alt = CONF["data_synonym"][name]
    if len(alt) < len(name):
        return alt
    return name


def to_data_name(name):
    alt = CONF["data_synonym"][name]
    if len(alt) < len(name):
        return name
    return alt


def get_short_data_names():
    return CONF["short_data_names"]
