from os.path import exists

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
    if _check_has_set_files(filepath, ["bed", "bim", "fam"]):
        return "bed"
    if any([filepath.endswith(ext) for ext in imexts]):
        return "image"
    if filepath.endswith(".txt"):
        return "csv"
    if filepath.endswith(".bgen"):
        return "bgen"
    if filepath.endswith(".gemma"):
        return "bimbam-pheno"
    if _check_has_set_files(filepath, ["rel", "rel.id"]):
        return "plink2-rel"
    return "unknown"


def _check_has_set_files(filepath, exts):
    files = [filepath + "." + ext for ext in exts]
    ok = [exists(f) for f in files]

    if sum(ok) > 0 and sum(ok) < len(exts):
        mfiles = ", ".join([files[i] for i in range(3) if not ok[i]])
        print("The following file(s) are missing:", mfiles)
        return False

    return all(ok)
