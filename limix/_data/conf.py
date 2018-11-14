CONF = {
    "dim_order": {
        "sample": 0,
        "trait": 1,
        "candidate": 1,
        "covariate": 1,
        "sample_0": 0,
        "sample_1": 1,
    },
    "dim_names": {"sample", "candidate", "covariate", "trait"},
    "data_names": {"trait", "genotype", "covariates", "covariance"},
    "short_data_names": {"y", "G", "M", "K"},
    "data_synonym": {
        "y": "trait",
        "trait": "y",
        "G": "genotype",
        "genotype": "G",
        "M": "covariates",
        "covariates": "M",
        "K": "covariance",
        "covariance": "K",
    },
    "data_dims": {"trait": ["sample", "trait"], "genotype": ["sample", "candidate"]},
}


def get_data_dims(name):
    return CONF["data_dims"][name]


def is_dim_name(name):
    return name in CONF["dim_names"]


def is_data_name(name):
    return name in CONF["data_names"]


def is_short_data_name(name):
    return name in CONF["short_data_names"]


def dim_order(name):
    return CONF["dim_order"][name]


def short_data_name(name):
    alt = CONF["data_synonym"][name]
    if len(alt) < len(name):
        return alt
    return name


def data_name(name):
    alt = CONF["data_synonym"][name]
    if len(alt) < len(name):
        return name
    return alt


def short_data_names():
    return CONF["short_data_names"]


def is_dim_hint(dim):
    return dim[0] == "_" and is_dim_name(dim[1:])


def dim_name_to_hint(dim):
    if not (is_dim_name(dim) or is_dim_hint(dim)):
        raise ValueError("`{}` is not a dimensional name/hint.")
    if is_dim_name(dim):
        return "_" + dim
    return dim


def dim_hint_to_name(dim):
    if not (is_dim_name(dim) or is_dim_hint(dim)):
        raise ValueError("`{}` is not a dimensional name/hint.")
    if is_dim_hint(dim):
        return dim[1:]
    return dim
