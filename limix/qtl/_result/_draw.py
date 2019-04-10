from ._table import Table


def draw_title(title):
    msg = f"{title}\n"
    msg += "=" * len(title) + "\n\n"
    return msg


def draw_model(lik, mean, covariance):
    if lik == "normal":
        var = "𝐲"
    else:
        var = "𝐳"

    msg = f"{var} ~ 𝓝({mean}, {covariance})"
    msg += _lik_formulae(lik)
    return msg


def draw_alt_hyp_table(hyp_num, stats, effsizes):
    cols = ["lml", "cov. effsizes", "cand. effsizes"]
    table = Table(cols, index=_describe_index())
    table.add_column(_describe(stats, f"lml{hyp_num}"))
    df = effsizes[f"h{hyp_num}"]
    table.add_column(_describe(df[df["effect_type"] == "covariate"], "effsize"))
    table.add_column(_describe(df[df["effect_type"] == "candidate"], "effsize"))
    return "\n" + table.draw() + "\n\n"


def draw_lrt_table(test_titles, pv_names, stats):
    table = Table(test_titles, index=_describe_index())

    for name in pv_names:
        pv = stats[name].describe().iloc[1:]
        table.add_column(pv)

    return table.draw()


def _lik_formulae(lik):
    msg = ""

    if lik == "bernoulli":
        msg += f" for yᵢ ~ Bern(μᵢ=g(zᵢ)) and g(x)=1/(1+e⁻ˣ)\n"
    elif lik == "probit":
        msg += f" for yᵢ ~ Bern(μᵢ=g(zᵢ)) and g(x)=Φ(x)\n"
    elif lik == "binomial":
        msg += f" for yᵢ ~ Binom(μᵢ=g(zᵢ), nᵢ) and g(x)=1/(1+e⁻ˣ)\n"
    elif lik == "poisson":
        msg += f" for yᵢ ~ Poisson(λᵢ=g(zᵢ)) and g(x)=eˣ\n"
    else:
        msg += "\n"
    return msg


def _describe(df, field):
    return df[field].describe().iloc[1:]


def _describe_index():
    return ["mean", "std", "min", "25%", "50%", "75%", "max"]
