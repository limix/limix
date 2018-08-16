# import numpy as np


# def _impl(a, func_name):
#     import dask.array as da

#     if isinstance(a, da.Array):
#         return getattr(da, func_name)
#     return getattr(np, func_name)


# def log(a, *args, **kwargs):
#     f = _impl(a, "log")
#     return f(a, *args, **kwargs)


# def isnan(a, *args, **kwargs):
#     f = _impl(a, "isnan")
#     return f(a, *args, **kwargs)


# def zeros_like(a, *args, **kwargs):
#     f = _impl(a, "zeros_like")
#     return f(a, *args, **kwargs)


# def abs(a, *args, **kwargs):
#     f = _impl(a, "abs")
#     return f(a, *args, **kwargs)


# def nanmean(a, *args, **kwargs):
#     f = _impl(a, "nanmean")
#     return f(a, *args, **kwargs)


# def nanstd(a, *args, **kwargs):
#     f = _impl(a, "nanstd")
#     return f(a, *args, **kwargs)


# def clip(a, *args, **kwargs):
#     f = _impl(a, "clip")
#     return f(a, *args, **kwargs)


# def copyto(a, *args, **kwargs):
#     f = _impl(a, "copyto")
#     return f(a, *args, **kwargs)


# def zeros(a, *args, **kwargs):
#     f = _impl(a, "zeros")
#     return f(a, *args, **kwargs)


# def sum(a, *args, **kwargs):
#     f = _impl(a, "sum")
#     return f(a, *args, **kwargs)


# def cumsum(a, *args, **kwargs):
#     f = _impl(a, "cumsum")
#     return f(a, *args, **kwargs)


# def flipud(a, *args, **kwargs):
#     f = _impl(a, "flipud")
#     return f(a, *args, **kwargs)


# def log10(a, *args, **kwargs):
#     f = _impl(a, "log10")
#     return f(a, *args, **kwargs)


# def all(a, *args, **kwargs):
#     f = _impl(a, "all")
#     return f(a, *args, **kwargs)


# def isfinite(a, *args, **kwargs):
#     f = _impl(a, "isfinite")
#     return f(a, *args, **kwargs)


# def any(a, *args, **kwargs):
#     f = _impl(a, "any")
#     return f(a, *args, **kwargs)


# def asarray(a, *args, **kwargs):
#     import dask.array as da

#     if isinstance(a, da.Array):
#         return da.asarray(a)

#     return np.asarray(a, *args, **kwargs)
