def fetch(spec, target, verbose):
    import limix

    spec = _parse_fetch_spec(spec)
    if verbose:
        utarget = target[0].upper() + target[1:]
        print("{} file type: {}".format(utarget, spec["filetype"]))

    return limix.io.fetch(target, spec, verbose=verbose)


def _parse_fetch_spec(fetch_spec):
    from ..io._detect import infer_filetype

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
