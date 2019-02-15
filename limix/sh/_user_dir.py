def user_cache_dir():
    import os
    from appdirs import user_cache_dir as ucdir

    d = ucdir("limix")
    os.makedirs(d, exist_ok=True)
    return d
