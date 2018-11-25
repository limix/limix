CONF = {
    "dim_axis": {
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
