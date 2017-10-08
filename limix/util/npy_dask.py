# TODO: document this module
def unary_choice(unary_func_name):

    def func(a, *args, **kwargs):
        import dask.array as da
        import numpy as np

        if isinstance(a, da.Array):
            f = getattr(da, unary_func_name)
        else:
            f = getattr(np, unary_func_name)

        return f(a, *args, **kwargs)

    return func


funcs = ['log', 'isnan', 'zeros_like', 'abs', 'nanmean', 'nanstd',
         'clip', 'nanmean', 'copyto', 'zeros', 'sum',
         'cumsum', 'flipud', 'log10', 'all', 'isfinite', 'any']

g = globals()
for f in funcs:
    g[f] = unary_choice(f)


def asarray(a, *args, **kwargs):
    import dask.array as da
    import numpy as np

    if isinstance(a, da.Array):
        return da.asarray(a)

    return np.asarray(a, *args, **kwargs)
