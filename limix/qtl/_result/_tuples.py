from collections import namedtuple

VariantResult = namedtuple(
    "VariantResult",
    [
        "lml",
        "covariate_effsizes",
        "candidate_effsizes",
        "covariate_effsizes_se",
        "candidate_effsizes_se",
        "scale",
    ],
)

Result = namedtuple("Result", ["idx", "h1", "h2"])
