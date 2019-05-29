from click import UsageError
import os


class InputData:
    def __init__(self):
        self._types = {
            "phenotype-matrices": [],
            "covariates-matrix": None,
            "genotype-matrix": None,
            "kinship-matrix": None,
        }

    def process_plink1(self, params):
        exts = set(["bim", "fam", "bed"])
        wanted = [x for x in params if x[0] in exts]

        if len([w[0] for w in wanted]) > len(set([w[0] for w in wanted])):
            msg = "The --bim, --fam, and --bed options cannot be repeated."
            raise UsageError(msg)

        if len(wanted) < len(exts):
            msg = "If you define one of the --bim, --fam, and --bed options, you"
            msg += " must also define the other two."
            raise UsageError(msg)

        wanted = dict(wanted)
        self._read_plink1_bin(wanted["bed"], wanted["bim"], wanted["fam"])

        remain = [x for x in params if x[0] not in exts]

        return remain

    def _read_plink1_bin(self, bed, bim, fam):
        from pandas_plink import read_plink1_bin

        G = read_plink1_bin(bed, bim, fam, verbose=False)

        if self._types["genotype-matrix"] is not None:
            raise UsageError("More than one genotype matrix has been defined.")
        self._types["genotype-matrix"] = G

        y = G["trait"]
        for field in [f for f in y.coords.keys() if f != y.dims[0]]:
            del y[field]

        y.name = os.path.basename(fam) + "_trait"
        self._types["phenotype-matrices"].append(y)

    def process_grm(self, filespec):
        from xarray import DataArray
        from numpy import loadtxt

        filespec = filespec.split(":")

        if len(filespec) == 1:
            filespec.append("unknown")

        if len(filespec) > 2:
            raise ValueError()

        K = DataArray(loadtxt(filespec[0]), dims=["sample_0", "sample_1"])
        self._types["kinship-matrix"] = K

    @property
    def phenotypes(self):
        return self._types["phenotype-matrices"]

    @property
    def covariates(self):
        return self._types["covariates-matrix"]

    @property
    def genotype(self):
        return self._types["genotype-matrix"]

    @property
    def kinship(self):
        return self._types["kinship-matrix"]


class ProcessInputData:
    def __init__(self, params):
        self._params = params
        self._input_data = InputData()
        self._look_for_kinship()
        self._look_for_genotype()
        self._look_for_covariates()
        self._look_for_phenotypes()

    def _look_for_kinship(self):
        remain = []
        params = self._params

        while len(params) > 0:
            if params[0][0] == "grm":
                self._input_data.process_grm(params[0][1])
            else:
                remain.append(params[0])
            del params[0]

        self._params = remain

    def _look_for_genotype(self):
        remain = []
        params = self._params

        while len(params) > 0:
            if _is_plink1_bin_option(params[0]):
                params = self._input_data.process_plink1(params)
            else:
                remain.append(params[0])
                del params[0]

        self._params = remain

    def _look_for_covariates(self):
        pass

    def _look_for_phenotypes(self):
        pass

    @property
    def input_data(self):
        return self._input_data


def _is_plink1_bin_option(param):
    return param[0] in set(["bim", "fam", "bed"])


def _bed_option(filepath, **_):
    from limix.io.plink import _read_bed

    return _read_bed(filepath)
