def is_array(a):
    pkg = a.__class__.__module__.split(".")[0]
    name = a.__class__.__name__

    return pkg == "numpy" and name == "ndarray"
