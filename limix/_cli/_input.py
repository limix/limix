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
        self._pheno_names = None

    def set_opt(self, opt, **kwargs):
        {
            "pheno": self._set_pheno,
            "pheno-name": self._set_pheno_name,
            "bfile": self._set_bfile,
            "bed": self._set_bed,
            "grm": self._set_grm,
            "rel": self._set_rel,
        }[opt](**kwargs)

    def _set_pheno(self, filepath):
        from limix.io import plink

        if filepath is None:
            return

        y = plink.read_pheno(filepath)
        y.name = _trait_name_from_filepath(filepath)

        self._types["phenotype-matrices"].append(y)

    def _set_pheno_name(self, pheno_name):
        if pheno_name is None:
            self._pheno_names = None
        else:
            self._pheno_names = [n.strip() for n in pheno_name.split(",")]

    def _set_bfile(self, bfile_prefix):
        from pandas_plink import read_plink1_bin

        if bfile_prefix is None:
            return

        self._types["genotype-matrix"] = read_plink1_bin(bfile_prefix + ".bed")

    def _set_bed(self, bed_filepath, bim_filepath, fam_filepath):
        from pandas_plink import read_plink1_bin

        if bed_filepath is None:
            return

        G = read_plink1_bin(bed_filepath, bim_filepath, fam_filepath, verbose=False)

        if self._types["genotype-matrix"] is not None:
            raise UsageError("More than one genotype matrix has been defined.")
        self._types["genotype-matrix"] = G

        y = G["trait"]
        for field in [f for f in y.coords.keys() if f != y.dims[0]]:
            del y[field]

        y.name = _trait_name_from_filepath(fam_filepath)
        self._types["phenotype-matrices"].append(y)

    def _set_grm(self, filepath):
        from pandas_plink import read_grm

        if filepath is None:
            return

        self._types["kinship-matrix"] = read_grm(filepath)

    def _set_rel(self, filepath):
        from pandas_plink import read_rel

        if filepath is None:
            return

        self._types["kinship-matrix"] = read_rel(filepath)

    @property
    def phenotypes(self):
        import click
        from xarray import concat

        if len(self._types["phenotype-matrices"]) == 0:
            return None

        Y = concat(self._types["phenotype-matrices"], dim="trait")
        if self._pheno_names is not None:
            try:
                Y = Y.sel(trait=self._pheno_names)
            except KeyError:
                read_traits = Y.trait.values.tolist()
                msg = "not all specified phenotypes have been found.\n"
                msg += f"Specified phenotypes: {self._pheno_names}\n"
                msg += f"Loaded phenotypes: {read_traits}\n"
                raise click.UsageError(msg)

        return Y

    @property
    def covariates(self):
        return self._types["covariates-matrix"]

    @property
    def genotype(self):
        return self._types["genotype-matrix"]

    @property
    def kinship(self):
        return self._types["kinship-matrix"]

    def __str__(self):
        from limix._display import summarize_list_repr

        msg = "Requested phenotypes: "
        if self._pheno_names is None:
            msg += "<all phenotypes>"
        else:
            msg += summarize_list_repr(self._pheno_names, 5)
        msg += "\n\n"

        msg += _repr_input("Phenotype", self.phenotypes) + "\n" * 2
        msg += _repr_input("Covariate", self.covariates) + "\n" * 2
        msg += _repr_input("Genotype", self.genotype) + "\n" * 2
        msg += _repr_input("Kinship", self.kinship)
        return msg


def _repr_input(header, X):
    msg = f"{header}\n"
    msg += ("-" * len(header)) + "\n"
    if X is None:
        name = header.lower()
        msg += f"No {name} has been provided."
    else:
        msg += str(X)
    return msg


class ProcessInputData:
    def __init__(self, params):
        self._params = params
        self._input_data = QTLInputData()
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
