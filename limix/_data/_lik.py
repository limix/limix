def check_likelihood_name(likname):

    valid_names = set(["normal", "bernoulli", "probit", "binomial", "poisson"])
    if likname not in valid_names:
        msg = "Unrecognized likelihood name: {}.\n".format(likname)
        msg += "Valid names are: {}.".format(valid_names)
        raise ValueError(msg)
