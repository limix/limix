def is_series(a):
    pkg = a.__class__.__module__.split(".")[0]
    name = a.__class__.__name__

    return pkg == "pandas" and name == "Series"


def is_data_frame(a):
    pkg = a.__class__.__module__.split(".")[0]
    name = a.__class__.__name__

    return pkg == "pandas" and name == "DataFrame"
