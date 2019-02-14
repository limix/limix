from collections import namedtuple

FetchSpec = namedtuple("FetchSpec", ["filepath", "filetype", "matrix_spec"])


def fetch(target, fetch_spec, verbose=True):
    from .._data import asarray
    from .._data import assert_target, assert_filetype

    assert_target(target)

    if isinstance(fetch_spec, str):
        fetch_spec = _parse_fetch_spec(fetch_spec)

    assert_filetype(fetch_spec.filetype)
    matrix_spec = fetch_spec.matrix_spec
    dims = {d: matrix_spec[d] for d in ["row", "col"] if d in matrix_spec}

    X = _dispatch[target][fetch_spec.filetype](fetch_spec.filepath, verbose=verbose)
    X = asarray(X, target, dims)

    if len(matrix_spec["sel"]) > 0:
        X = X.sel(**matrix_spec["sel"])

    X = asarray(X, target)
    X.name = target

    return X


def _fetch_npy_covariance(filepath, verbose=True):
    from .npy import read

    return read(filepath, verbose=verbose)


def _fetch_bed_genotype(filepath, verbose=True):
    from .plink import read
    from limix._data import asarray

    candidates, samples, G = read(filepath, verbose=verbose)

    assert G.shape == (len(candidates), len(samples))
    G = asarray(G, "genotype", ["candidate", "sample"])
    # Make sure we return the same matrix layout that has been read.
    G = G.transpose("candidate", "sample")

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


_dispatch = {
    "genotype": {"bed": _fetch_bed_genotype},
    "trait": {"bimbam-pheno": _fetch_bimbam_phenotype, "csv": _fetch_csv_phenotype},
    "covariance": {"npy": _fetch_npy_covariance},
}


def _parse_fetch_spec(spec):
    import os
    from ._detect import infer_filetype

    drive, spec = os.path.splitdrive(spec)

    ncolons = sum([1 for c in spec if c == ":"])
    if ncolons > 2:
        raise ValueError("Wrong specification syntax: there were more than two colons.")

    spec = drive + spec
    spec = spec + ":" * (2 - ncolons)
    spec = dict(zip(["filepath", "filetype", "matrix_spec"], spec.split(":")))

    if spec["filetype"] == "":
        spec["filetype"] = infer_filetype(spec["filepath"])
    spec["matrix_spec"] = _parse_matrix_spec(spec["matrix_spec"])

    return FetchSpec(**spec)


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
