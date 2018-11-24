from functools import wraps


def return_none_if_none(f):
    """ Return ``None`` if first argument is ``None``. """

    @wraps(f)
    def new_func(x, *args, **kwargs):
        if x is None:
            return None
        return f(x, *args, **kwargs)

    return new_func
