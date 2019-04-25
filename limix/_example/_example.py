_filenames = {
    "data.csv": {
        "url": "http://rest.s3for.me/limix/example/data.csv",
        "sha256": "9b09e5d569b256aabc42d3e04139f7f55789283ee642023adee455b6e5ef7d14",
    },
    "data.h5": {
        "url": "http://rest.s3for.me/limix/example/data.h5",
        "sha256": "30bc0374211f4ad7f763f631fc7cf8d1fdee78d5a0064facc8fcf09ab9ace64c",
    },
    "example.gen": {
        "url": "http://rest.s3for.me/limix/example/example.gen",
        "sha256": "a21903b261254981edf552798bda8f911423bf92083759b93f3e78295cc7e382",
    },
    "example.sample": {
        "url": "http://rest.s3for.me/limix/example/example.sample",
        "sha256": "f10f3304104fbecfd15703da06ef3e9d922a5dc069644bacb3f952d568d4c9db",
    },
    "expr.csv": {
        "url": "http://rest.s3for.me/limix/example/expr.csv",
        "sha256": "48ab488b1a1a7651d1152dae1e13605df9a342667ccba48328ae4b86b3ac359a",
    },
    "pheno.csv": {
        "url": "http://rest.s3for.me/limix/example/pheno.csv",
        "sha256": "9b09e5d569b256aabc42d3e04139f7f55789283ee642023adee455b6e5ef7d14",
    },
    "xarr.hdf5": {
        "url": "http://rest.s3for.me/limix/example/xarr.hdf5",
        "sha256": "fde8e4830d04d47f4463ba17704679b3259dc4907f56cc2d69a930a584ea6a00",
    },
    "phenotype.gemma": {
        "url": "http://rest.s3for.me/limix/example/phenotype.gemma",
        "sha256": "90089e70364dc5e2b6d71ea00e3b5e60203d9940ca909cf2ce688aca5178044b",
    },
    "data.bed": {
        "url": "http://rest.s3for.me/limix/example/data.bed",
        "sha256": "7bf9162d486116db12686e14ef7bf2d53c46bb788041555041891cebabf1ba09",
    },
    "data.bim": {
        "url": "http://rest.s3for.me/limix/example/data.bim",
        "sha256": "f6b70bebe81aad8941903cd2a53bcf30507efdd7c52d6692e8d62cb90a7ddb1c",
    },
    "data.fam": {
        "url": "http://rest.s3for.me/limix/example/data.fam",
        "sha256": "701a94fc09860bc341b42497ac9cae893e771798270d69d3fb566cfa09c03c8f",
    },
    "data.h5.bz2": {
        "url": "http://rest.s3for.me/limix/example/data.h5.bz2",
        "sha256": "2df755158b403416e1b0270815b2fc18f7d23ec5590c9cabacd2e1be98510c59",
    },
    "chrom22_subsample20_maf0.10.bed": {
        "url": "http://rest.s3for.me/limix/example/chrom22_subsample20_maf0.10.bed",
        "sha256": "cf18099968dbb9a514d2e9b86f02f2fab2c549e7859d7560ae5bc79a841dd269",
    },
    "chrom22_subsample20_maf0.10.fam": {
        "url": "http://rest.s3for.me/limix/example/chrom22_subsample20_maf0.10.fam",
        "sha256": "76bd241625bf2b6aa5fe763d4d037609fbd8c70eadf261faeab99cdd5fe7443e",
    },
    "chrom22_subsample20_maf0.10.bim": {
        "url": "http://rest.s3for.me/limix/example/chrom22_subsample20_maf0.10.bim",
        "sha256": "24824433bef33e17a9f9713868f492c9957efbc585e0010dfcc8536fdf7996d4",
    },
    "expr_nan.csv": {
        "url": "http://rest.s3for.me/limix/example/expr_nan.csv",
        "sha256": "6ec7d2a21dbb2b2e65b0a28cf7f7bdbf6cbb974f9cfd1d653050e6b1aa84dfbb",
    },
}


def _download(url, filepath, hash):
    import os
    from urllib.request import urlretrieve
    import limix

    local_filename = urlretrieve(url)[0]
    try:
        os.rename(local_filename, filepath)
    except OSError:
        import shutil

        shutil.move(local_filename, filepath)

    if limix.sh.filehash(filepath) != hash:
        raise RuntimeError("Hash does not match for file {}.".format(filepath))


class file_example(object):
    """
    Put file examples into a temporary folder.
    """

    def __init__(self, filenames):
        import tempfile
        from limix.sh._dir import makedirs
        import os
        from os.path import join
        from limix.sh._user_dir import user_cache_dir
        import limix

        self._unlist = False
        if not isinstance(filenames, (tuple, list)):
            filenames = [filenames]
            self._unlist = True

        self._orig_folder = join(user_cache_dir(), "examples")
        makedirs(self._orig_folder)

        for fn in filenames:
            if fn not in _filenames:
                raise ValueError(f"Unrecognized file name {fn}.")
            info = _filenames[fn]
            fp = join(self._orig_folder, fn)
            if os.path.exists(fp):
                if limix.sh.filehash(fp) != info["sha256"]:
                    os.remove(fp)
                    _download(info["url"], fp, info["sha256"])
            else:
                _download(info["url"], fp, info["sha256"])

        self._dirpath = tempfile.mkdtemp()
        self._filenames = filenames

    def __enter__(self):
        import shutil
        from os.path import join

        filepaths = []
        for fn in self._filenames:
            fp = join(self._dirpath, fn)
            shutil.copy(join(self._orig_folder, fn), fp)
            filepaths.append(fp)

        if self._unlist:
            return filepaths[0]
        return filepaths

    def __exit__(self, *_):
        import shutil

        try:
            shutil.rmtree(self._dirpath)
        except PermissionError:
            pass
