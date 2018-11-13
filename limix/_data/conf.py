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
}


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
