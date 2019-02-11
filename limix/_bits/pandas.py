__all__ = ["is_series", "is_dataframe"]


def is_series(a):
    pkg = a.__class__.__module__.split(".")[0]
    name = a.__class__.__name__

    return pkg == "pandas" and name == "Series"


def is_dataframe(a):
    pkg = a.__class__.__module__.split(".")[0]
    name = a.__class__.__name__

    return pkg == "pandas" and name == "DataFrame"
