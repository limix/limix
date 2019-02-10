from ._pipeline import _Pipeline


def _preprocessing(
    data, filter_specs, filter_missing, filter_maf, impute, normalize, verbose
):
    # from limix._data import conform_dataset
    # from limix._display import session_line

    pipeline = _Pipeline(data)

    for spec in filter_specs:
        # data = _process_filter(f, data)
        for target in data.keys():
            pipeline.append(_process_filter, spec, target)
        # for target in data.keys():
        #     layout.append(target, "filter {}".format(i), data[target].shape)
        #     if data["y"].sample.size == 0:
        #         print(layout.to_string())
        #         raise RuntimeError("Exiting early because there is no sample left.")

    # for i, f in enumerate(filter):
    #     data = _process_filter(f, data)
    #     for target in data.keys():
    #         layout.append(target, "filter {}".format(i), data[target].shape)
    #         if data["y"].sample.size == 0:
    #             print(layout.to_string())
    #             raise RuntimeError("Exiting early because there is no sample left.")

    # for f in filter_missing:
    #     with session_line("Applying `{}`... ".format(f)):
    #         _process_filter_missing(f, data)
    #         if data["y"].sample.size == 0:
    #             print(layout.to_string())
    #             raise RuntimeError("Exiting early because there is no sample left.")

    # if filter_maf is not None:
    #     with session_line("Removing candidates with MAF<{}... ".format(filter_maf)):
    #         data["G"] = _process_filter_maf(float(filter_maf), data["G"])

    #     for target in data.keys():
    #         layout.append(target, "maf filter", data[target].shape)

    #     if data["G"].candidate.size == 0:
    #         print(layout.to_string())
    #         raise RuntimeError("Exiting early because there is no candidate left.")

    # for imp in impute:
    #     with session_line("Imputting missing values (`{}`)... ".format(imp)):
    #         data = _process_impute(imp, data)

    for spec in normalize:
        for target in data.keys():
            pipeline.append(_process_normalize, spec, target)

    pipeline.run()

    return data


def _process_filter(data, layout, filter_spec, target):
    from limix._bits.xarray import query
    from limix._data import to_short_data_name

    elems = [fs.strip() for fs in filter_spec.strip().split(":")]
    if len(elems) != 2:
        msg = "Filter syntax error. It should have been\n"
        msg += "  <TARGET>: <COND>\n"
        msg += f"but we received `{filter_spec}`."
        raise ValueError(msg)

    target_name = elems[0]
    expr = elems[1]

    if to_short_data_name(target_name) == target:
        data[target] = query(data[target], expr)

    layout.append(target, "filter", data[target].shape)

    return data


def _process_normalize(data, layout, normalize_spec, target_short_name):
    import limix
    from limix._data import to_short_data_name, dim_hint_to_name, dim_name_to_hint

    elems = [fs.strip() for fs in normalize_spec.strip().split(":")]
    if len(elems) != 3:
        msg = "Normalize syntax error. It should have been\n"
        msg += "  <TARGET>:<DIM>:<METHOD>\n"
        msg += f"but we received `{normalize_spec}`."
        raise ValueError(msg)

    target = elems[0]
    dim = elems[1]
    method = elems[2]

    target = to_short_data_name(elems[0])
    if len(elems[1]) == 0:
        dim = elems[0]
    else:
        dim = elems[1]

    if len(elems[2]) > 0:
        method = elems[2]
    else:
        method = "gaussianize"

    def in_dim(X, dim):
        return dim_hint_to_name(dim) in X.dims or dim_name_to_hint(dim) in X.dims

    X = data[target]
    if not in_dim(X, dim):
        raise ValueError("Unrecognized dimension: {}.".format(dim))

    axis = next(i for i in range(len(X.dims)) if in_dim(X, dim))

    if method == "gaussianize":
        if axis == 0:
            X = limix.qc.quantile_gaussianize(X.T).T
        else:
            X = limix.qc.impute.mean_impute(X)
    else:
        raise ValueError("Unrecognized normalization method: {}.".format(method))

    data[target] = X

    return data


def _process_filter_missing(expr, data):
    elems = [e.strip() for e in expr.strip().split(":")]
    if len(elems) < 2 or len(elems) > 3:
        raise ValueError("Missing filter syntax error.")

    target = elems[0]
    dim = elems[1]

    if len(elems) == 3:
        how = elems[2]
    else:
        how = "any"

    data[target] = data[target].dropna(dim, how)


def _process_filter_maf(maf, G):
    from limix import compute_maf

    mafs = compute_maf(G)
    ok = mafs >= maf
    return G.isel(candidate=ok)


def _process_impute(expr, data):
    import limix
    from limix._data import to_short_data_name, dim_hint_to_name, dim_name_to_hint

    elems = [e.strip() for e in expr.strip().split(":")]
    if len(elems) < 2 or len(elems) > 3:
        raise ValueError("Imputation syntax error.")

    target = to_short_data_name(elems[0])
    dim = elems[1]

    if len(elems) == 3:
        method = elems[2]
    else:
        method = "mean"

    def in_dim(X, dim):
        return dim_hint_to_name(dim) in X.dims or dim_name_to_hint(dim) in X.dims

    X = data[target]
    if not in_dim(X, dim):
        raise ValueError("Unrecognized dimension: {}.".format(dim))

    if method == "mean":
        axis = next(i for i in range(len(X.dims)) if in_dim(X, dim))
        if axis == 0:
            X = limix.qc.impute.mean_impute(X.T).T
        else:
            X = limix.qc.impute.mean_impute(X)
    else:
        raise ValueError("Unrecognized imputation method: {}.".format(method))

    data[target] = X

    return data
