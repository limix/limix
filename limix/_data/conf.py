CONF = {
    "dim_order": {"sample": 0, "trait": 1, "candidate": 1, "covariate": 1},
    "dim_names": {"sample", "candidate", "covariate", "trait"},
    "data_names": {"trait", "genotype", "covariates", "covariance"},
    "data_synonym": {"y": "trait", "trait": "y", "G": "genotype", "genotype": "G"},
}


def is_dim_name(name):
    return name in CONF["dim_names"]


def is_data_name(name):
    return name in CONF["data_names"]


def dim_order(name):
    return CONF["dim_order"][name]


def short_data_name(name):
    alt = CONF["synonym"][name]
    if len(alt) < len(name):
        return alt
    return name


def long_data_name(name):
    alt = CONF["synonym"][name]
    if len(alt) < len(name):
        return name
    return alt
