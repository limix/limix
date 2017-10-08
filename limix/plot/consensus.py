from __future__ import absolute_import, division

from numpy import asarray, logical_and, ones, sqrt


# TODO: finish the documentation
class ConsensusCurve(object):
    r"""Consolidate multiple curves in a single one."""

    def __init__(self):
        self._x = []
        self._y = []

    def add(self, x, y):
        r"""Add a new curve."""
        self._x.append(asarray(x))
        self._y.append(asarray(y))

    def consensus(self, std_dev=3.0):
        r"""Return a consensus curve."""
        x = self._x
        y = self._y
        x, ybottom, y, ytop = _consensus_curve(x, y, std_dev=std_dev)

        return (x, ybottom, y, ytop)


def _2dbulk(x, left, right, y):
    xbulk = []
    ybulk = []
    for (i, xi) in enumerate(x):
        ok = logical_and(left <= xi, xi <= right)
        xbulk.append(xi[ok].mean())
        ybulk.append(y[i][ok].mean())

    return (asarray(xbulk), asarray(ybulk))


def _consensus_curve(x, y, std_dev=3.0):
    n = -1
    for xi in x:
        n = max(len(xi), n)

    nx = ones((len(x), n))
    ny = ones((len(y), n))

    for i, xi in enumerate(x):
        nx[i, :len(xi)] = xi

    for i, yi in enumerate(y):
        ny[i, :len(yi)] = yi

    x = nx
    y = ny

    xavg = x.mean(0)
    yavg = y.mean(0)

    err = y.std(0) * std_dev / sqrt(x.shape[0])

    return (xavg, yavg - err, yavg, yavg + err)
