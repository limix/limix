def check_likelihood_name(likname):
    from ._conf import CONF

    if likname not in CONF["likelihoods"]:
        msg = "Unrecognized likelihood name: {}.\n".format(likname)
        msg += "Valid names are: {}.".format(list(CONF["likelihoods"]))
        raise ValueError(msg)
