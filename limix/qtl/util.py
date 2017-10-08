from __future__ import division


def print_analysis(lik, name):
    lik_name = lik[0].upper() + lik[1:]
    print("*** %s using %s-GLMM ***" % (name, lik_name))
