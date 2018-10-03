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


def detect_file_type(filepath):
    # TODO document

    spec = _get_file_type_spec(filepath)
    if spec is not None:
        if spec in recognized_file_types:
            return spec

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
    return "unknown"


def _is_bed(filepath):
    files = [filepath + ext for ext in [".bed", ".bim", ".fam"]]
    ok = [exists(f) for f in files]

    if sum(ok) > 0 and sum(ok) < 3:
        mfiles = ", ".join([files[i] for i in range(3) if not ok[i]])
        print("The following file(s) are missing:", mfiles)
        return False

    return all(ok)


def _get_file_type_spec(filepath):
    filename = basename(filepath)
    if ":" not in filename:
        return None

    return filename.split(":")[-1]
