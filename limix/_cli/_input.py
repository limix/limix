import os

from click import UsageError


class QTLInputData:
    """
    Store user-provided, command-line input for QTL analysis.
    """

    def __init__(self):
        self._types = {
            "trait-matrices": [],
            "covariates-matrix": None,
            "genotype-matrix": None,
            "kinship-matrix": None,
        }
        self._trait_names = None

    def set_opt(self, opt, **kwargs):
        {
            "trait": self._set_trait,
            "trait-name": self._set_trait_name,
            "bfile": self._set_bfile,
            "bed": self._set_bed,
            "grm": self._set_grm,
            "rel": self._set_rel,
            "outdir": self._set_outdir,
        }[opt](**kwargs)

    def _set_trait(self, filepath):
        from limix.io import csv
        import xarray as xr

        if filepath is None:
            return

        y = csv.read(filepath)
        y = y.set_index("sample", drop=True)
        y = xr.DataArray(y, dims=["sample", "trait"])

        self._types["trait-matrices"].append(y)

    def _set_trait_name(self, trait_name):
        if trait_name is None:
            self._trait_names = None
        else:
            self._trait_names = [n.strip() for n in trait_name.split(",")]

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
        self._types["trait-matrices"].append(y)

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

    def _set_outdir(self, outdir):
        from os import makedirs

        if not outdir.exists():
            makedirs(outdir, exist_ok=True)

        self._outdir = outdir

    @property
    def traits(self):
        import click
        from xarray import concat

        if len(self._types["trait-matrices"]) == 0:
            return None

        Y = concat(self._types["trait-matrices"], dim="trait")
        if self._trait_names is not None:
            try:
                Y = Y.sel(trait=self._trait_names)
            except KeyError:
                read_traits = Y.trait.values.tolist()
                msg = "not all specified traits have been found.\n"
                msg += f"Specified traits: {self._trait_names}\n"
                msg += f"Loaded traits: {read_traits}\n"
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

    @property
    def outdir(self):
        return self._outdir

    def __str__(self):
        from limix._display import summarize_list_repr

        msg = "Requested traits: "
        if self._trait_names is None:
            msg += "<all traits>"
        else:
            msg += summarize_list_repr(self._trait_names, 5)
        msg += "\n\n"

        msg += _repr_input("Trait", self.traits) + "\n" * 2
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


def _trait_name_from_filepath(filepath):
    return os.path.splitext(os.path.basename(filepath))[0] + "_trait"
