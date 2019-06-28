import os

from click import UsageError


class QTLInputData:
    """
    Store user-provided, command-line input for QTL analysis.
    """

    def __init__(self):
        self._types = {
            "phenotype-matrices": [],
            "covariates-matrix": None,
            "genotype-matrix": None,
            "kinship-matrix": None,
        }

    def set_opt(self, opt, **kwargs):
        {
            "pheno": self.set_pheno,
            "bfile": self.set_bfile,
            "bed": self._set_bed,
            "grm": self.set_grm,
            "rel": self.set_rel,
        }[opt](**kwargs)

    def set_pheno(self, filepath):
        from limix.io import plink

        y = plink.read_pheno(filepath)
        y.name = _trait_name_from_filepath(filepath)

        self._types["phenotype-matrices"].append(y)

    def set_bfile(self, filepath):
        from pandas_plink import read_plink1_bin

        self._types["genotype-matrix"] = read_plink1_bin(filepath + ".bed")

    def set_bed(self, bed, bim, fam):
        from pandas_plink import read_plink1_bin

        G = read_plink1_bin(bed, bim, fam, verbose=False)

        if self._types["genotype-matrix"] is not None:
            raise UsageError("More than one genotype matrix has been defined.")
        self._types["genotype-matrix"] = G

        y = G["trait"]
        for field in [f for f in y.coords.keys() if f != y.dims[0]]:
            del y[field]

        y.name = _trait_name_from_filepath(fam)
        self._types["phenotype-matrices"].append(y)

    def set_rel(self, filespec):
        from pandas_plink import read_rel

        self._types["kinship-matrix"] = read_rel(filespec)

    def set_grm(self, filepath):
        from pandas_plink import read_grm

        self._types["kinship-matrix"] = read_grm(filepath)

    @property
    def phenotypes(self):
        from xarray import concat

        if len(self._types["phenotype-matrices"]) == 0:
            return None

        return concat(self._types["phenotype-matrices"], dim="trait")

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
            elif params[0][0] == "rel":
                self._input_data.process_rel(params[0][1])
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

        remain = []
        params = self._params

        while len(params) > 0:
            if params[0][0] == "pheno":
                self._input_data.process_pheno(params[0][1])
            else:
                remain.append(params[0])
            del params[0]

        self._params = remain

    @property
    def input_data(self):
        return self._input_data


def _is_plink1_bin_option(param):
    return param[0] in set(["bim", "fam", "bed", "bfile"])


def _bed_option(filepath, **_):
    from limix.io.plink import _read_bed

    return _read_bed(filepath)


def _trait_name_from_filepath(filepath):
    return os.path.splitext(os.path.basename(filepath))[0] + "_trait"
