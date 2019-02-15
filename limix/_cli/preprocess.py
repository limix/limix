def impute(data, spec):
    import limix
    from limix._data import CONF

    ncolons = sum(1 for c in spec if c == ":")
    if ncolons > 2:
        _syntax_error_msg("Impute", "<TARGET>:<DIM>:<METHOD>", spec)

    spec = spec.strip()
    spec = spec + ":" * (2 - ncolons)
    target, dim, method = [e.strip() for e in spec.split(":")]

    varname = CONF["target_to_varname"][target]
    x = data[varname]

    if dim == "":
        axis = -1
    else:
        axis = next(i for i, d in enumerate(x.dims) if d == dim)

    if method == "":
        method = "mean"

    if method == "mean":
        x = limix.qc.mean_impute(x, axis=axis)
    else:
        raise ValueError(f"Unrecognized impute method: {method}.")

    data[varname] = x

    return data


def normalize(data, spec):
    import limix
    from limix._data import CONF

    ncolons = sum(1 for c in spec if c == ":")
    if ncolons > 2:
        _syntax_error_msg("Normalize", "<TARGET>:<DIM>:<METHOD>", spec)

    spec = spec.strip()
    spec = spec + ":" * (2 - ncolons)
    target, dim, method = [e.strip() for e in spec.split(":")]

    varname = CONF["target_to_varname"][target]
    x = data[varname]

    if dim == "":
        axis = -1
    else:
        axis = next(i for i, d in enumerate(x.dims) if d == dim)

    if method == "":
        method = "gaussianize"

    if method == "gaussianize":
        x = limix.qc.quantile_gaussianize(x, axis=axis)
    elif method == "mean_std":
        x = limix.qc.mean_standardize(x, axis=axis)
    else:
        raise ValueError(f"Unrecognized normalization method: {method}.")

    data[varname] = x

    return data


def where(data, spec):
    from limix._bits.xarray import query
    from limix._data import CONF

    ncolons = sum(1 for c in spec if c == ":")
    if ncolons > 1:
        _syntax_error_msg("Where", "<TARGET>:<COND>", spec)

    spec = spec.strip()
    spec = spec + ":" * (1 - ncolons)
    target, cond = [e.strip() for e in spec.split(":")]

    varname = CONF["target_to_varname"][target]

    data[varname] = query(data[varname], cond)

    return data


def drop_missing(data, spec):
    from limix._data import CONF

    ncolons = sum(1 for c in spec if c == ":")
    if ncolons > 2:
        _syntax_error_msg("Drop-missing", "<TARGET>:<DIM>:<METHOD>", spec)

    spec = spec.strip()
    spec = spec + ":" * (2 - ncolons)
    target, dim, method = [e.strip() for e in spec.split(":")]

    if dim == "":
        if target == "trait":
            dim = "trait"
        elif target == "genotype":
            dim = "candidate"
        elif target == "covariate":
            dim = "covariate"

    if method == "":
        method = "any"

    varname = CONF["target_to_varname"][target]
    x = data[varname]

    if dim == "":
        dim = x.dims[-1]

    data[varname] = data[varname].dropna(dim, method)

    return data


def drop_maf(data, spec):
    raise NotImplementedError()


# def _process_filter_maf(maf, G):
#     from limix import compute_maf

#     mafs = compute_maf(G)
#     ok = mafs >= maf
#     return G.isel(candidate=ok)


def _syntax_error_msg(what, syntax, spec):
    msg = f"{what} syntax error. It should have been\n"
    msg += f"  {syntax}\n"
    msg += f"but we received `{spec}`."
    return msg
