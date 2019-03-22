def cache(func):
    from functools import lru_cache

    return lru_cache(maxsize=None)(func)
