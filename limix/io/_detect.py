from os.path import exists, basename

recognized_file_types = [
    "image",
    "hdf5",
    "csv",
    "npy",
    "grm.raw",
    "bed",
    "bgen",
    "bimbam-pheno",
]


def infer_filetype(filepath):
    imexts = [".png", ".bmp", ".jpg", "jpeg"]
    if filepath.endswith(".hdf5") or filepath.endswith(".h5"):
        return "hdf5"
    if filepath.endswith(".csv"):
        return "csv"
    if filepath.endswith(".npy"):
        return "npy"
    if filepath.endswith(".grm.raw"):
        return "grm.raw"
    if _is_bed(filepath):
        return "bed"
    if any([filepath.endswith(ext) for ext in imexts]):
        return "image"
    if filepath.endswith(".txt"):
        return "csv"
    if filepath.endswith(".bgen"):
        return "bgen"
    if filepath.endswith(".gemma"):
        return "bimbam-pheno"
    return "unknown"


def detect_filetype(fetch_spec):
    spec = _split_fetch_spec(fetch_spec)
    if spec["filetype"] != "":
        if spec["filetype"] in recognized_file_types:
            return spec["filetype"]
    return infer_filetype(spec["filepath"])


def get_fetch_spec(fetch_spec):
    spec = _split_fetch_spec(fetch_spec)
    if spec["filetype"] == "":
        spec["filetype"] = infer_filetype(spec["filepath"])
    spec["matrix_spec"] = _parse_matrix_spec(spec["matrix_spec"])
    return spec


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
            data["sel"].update({match.group(1): match.group(2)})

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


def _is_bed(filepath):
    files = [filepath + ext for ext in [".bed", ".bim", ".fam"]]
    ok = [exists(f) for f in files]

    if sum(ok) > 0 and sum(ok) < 3:
        mfiles = ", ".join([files[i] for i in range(3) if not ok[i]])
        print("The following file(s) are missing:", mfiles)
        return False

    return all(ok)

