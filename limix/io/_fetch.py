# from .._data import asarray


def fetch(target, fetch_spec, verbose=True):
    from .._data import asarray

    # from .._bits.xarray import hint_aware_sel

    # if not is_data_name(target):
    #     raise ValueError("`{}` is not a valid data name.".format(target))

    if not isinstance(fetch_spec, dict):
        fetch_spec = _parse_fetch_spec(fetch_spec)

    filetype = fetch_spec["filetype"]

    spec = fetch_spec["matrix_spec"]
    dims = {d: spec[d] for d in ["row", "col"] if d in spec}

    X = _dispatch[target][filetype](fetch_spec["filepath"], verbose=verbose)
    X = asarray(X, target, dims)
    # X = to_dataarray(X)
    # X = _read_dims_into(X, dims)

    if len(spec["sel"]) > 0:
        # X = hint_aware_sel(X, **spec["sel"])
        X = X.sel(**spec["sel"])

    X = asarray(X, target)
    X.name = target

    # if target == "trait":
    #     X = _set_missing_dim(X, get_dims_from_data_name(target)[: X.ndim])
    #     # breakpoint()
    #     # X = _sort_dims(X, get_dims_order_from_data_name(target))

    return X


def _fetch_npy_covariance(filepath, verbose=True):
    from .npy import read

    return read(filepath, verbose=verbose)


def _fetch_bed_genotype(filepath, verbose=True):
    from .plink import read
    from xarray import DataArray

    candidates, samples, G = read(filepath, verbose=verbose)

    G = DataArray(G.T)
    # G = DataArray(G.T, dims=get_dims_from_data_name("genotype"))

    for colname in samples.columns:
        G.coords[colname] = ("sample", samples[colname].values)

    G.coords[samples.index.name] = ("sample", samples.index)

    for colname in candidates.columns:
        G.coords[colname] = ("candidate", candidates[colname].values)

    G.coords[candidates.index.name] = ("candidate", candidates.index)

    return G


def _fetch_bimbam_phenotype(filepath, verbose):
    from .bimbam import read_phenotype

    return read_phenotype(filepath, verbose=verbose)


def _fetch_csv_phenotype(filepath, verbose):
    from .csv import read

    return read(filepath, verbose=verbose)


def _read_dims_into(X, dims):
    rc = {"row": 0, "col": 1}
    # If mentioned dims are already in the datarray, just transpose it
    # if necessary.
    for axis_name, dim_name in dims.items():
        try:
            first = next(i for i in range(len(X.dims)) if X.dims[i] == dim_name)
        except StopIteration:
            continue
        if rc[axis_name] != first:
            X = X.T

    # Se dim names if they were not found already in the dataarray.
    for axis_name, dim_name in dims.items():
        X = X.rename({X.dims[rc[axis_name]]: dim_name})
    return X


_dispatch = {
    "genotype": {"bed": _fetch_bed_genotype},
    "trait": {"bimbam-pheno": _fetch_bimbam_phenotype, "csv": _fetch_csv_phenotype},
    "covariance": {"npy": _fetch_npy_covariance},
}


def _fix_trait_dims(y):
    if y.ndim == 1:
        if y.dims[0] != "sample":
            y = y.rename({y.dims[0]: "sample"})
    else:
        y = _set_missing_dim(y, ["sample", "trait"])
    return y


def _set_missing_dim(arr, dims):
    unk_dims = set(arr.dims) - set(dims)
    if len(unk_dims) > 1:
        raise ValueError("Too many unknown dimension names.")
    elif len(unk_dims) == 1:
        known_dims = set(dims) - set(arr.dims)
        if len(known_dims) != 1:
            raise ValueError("Can't figure out what is the missing dimension name.")
        arr = arr.rename({unk_dims.pop(): known_dims.pop()})
    return arr


def _sort_dims(arr, dim_order):
    dim_order = {k: v for k, v in dim_order.items() if k in arr.dims}
    axis2dim = {v: k for k, v in dim_order.items()}
    axes = sorted(axis2dim.keys())
    return arr.transpose(*[axis2dim[a] for a in axes])


def _split_fetch_spec(txt):
    parts = []
    j = 0
    for i in range(len(txt)):
        if len(parts) == 2:
            parts.append(txt[i:])
        if txt[i] == ":":
            if j == i:
                raise ValueError("Invalid fetch specification syntax.")
            parts.append(txt[j:i])
            j = i + 1

    if len(txt[j:]) > 0:
        parts.append(txt[j:])

    if len(parts) == 0:
        raise ValueError("Invalid fetch specification syntax.")

    data = {"filepath": "", "filetype": "", "matrix_spec": ""}
    data["filepath"] = parts[0]
    if len(parts) > 1:
        data["filetype"] = parts[1]
    if len(parts) > 2:
        data["matrix_spec"] = parts[2]
    return data


def _parse_fetch_spec(fetch_spec):
    from ._detect import infer_filetype

    spec = _split_fetch_spec(fetch_spec)
    if spec["filetype"] == "":
        spec["filetype"] = infer_filetype(spec["filepath"])
    spec["matrix_spec"] = _parse_matrix_spec(spec["matrix_spec"])
    return spec


def _number_or_string(val):
    if "." in val:
        try:
            val = float(val)
            return val
        except ValueError:
            pass

    try:
        val = int(val)
        return val
    except ValueError:
        pass

    enclosed = False
    if val.startswith("'") and val.endswith("'") and len(val) > 1:
        enclosed = True

    if val.startswith("'") and val.endswith("'") and len(val) > 1:
        enclosed = True

    if not enclosed:
        val = '"' + val + '"'

    return val


def _parse_matrix_spec(txt):
    import re

    parts = _split_matrix_spec(txt)
    data = {"sel": {}}
    for p in parts:
        p = p.strip()
        if p.startswith("row"):
            data["row"] = p.split("=")[1]
        elif p.startswith("col"):
            data["col"] = p.split("=")[1]
        else:
            match = re.match(r"(^[^\[]+)\[(.+)\]$", p)
            if match is None:
                raise ValueError("Invalid fetch specification syntax.")
            # TODO: replace eval for something safer
            v = _number_or_string(match.group(2))
            data["sel"].update({match.group(1): eval(v)})

    return data


def _split_matrix_spec(txt):

    brackets = 0
    parts = []
    j = 0
    for i in range(len(txt)):
        if txt[i] == "[":
            brackets += 1
        elif txt[i] == "]":
            brackets -= 1
        elif txt[i] == "," and brackets == 0:
            if j == i:
                raise ValueError("Invalid fetch specification syntax.")
            parts.append(txt[j:i])
            j = i + 1
        if brackets < 0:
            raise ValueError("Invalid fetch specification syntax.")

    if len(txt[j:]) > 0:
        parts.append(txt[j:])

    return parts


# def hint_aware_sel(x, **kwargs):
#     from .._data import is_dim_hint, is_dim_name, dim_name_to_hint, dim_hint_to_name

#     for k in kwargs.keys():
#         if in_coords_dim(x, k):
#             continue
#         if is_dim_name(k) or is_dim_hint(k):
#             if in_coords_dim(x, dim_name_to_hint(k)):
#                 new_k = dim_name_to_hint(k)
#                 if new_k not in kwargs:
#                     kwargs[new_k] = kwargs[k]
#                     del kwargs[k]
#             elif in_coords_dim(x, dim_hint_to_name(k)):
#                 new_k = dim_hint_to_name(k)
#                 if new_k not in kwargs:
#                     kwargs[new_k] = kwargs[k]
#                     del kwargs[k]

#     return x.sel(**kwargs)
