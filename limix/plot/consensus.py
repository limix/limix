from __future__ import absolute_import, division

from numpy import asarray, logical_and, sqrt


class ConsensusCurve(object):
    def __init__(self):
        self._x = []
        self._y = []

    def add(self, x, y):
        self._x.append(asarray(x))
        self._y.append(asarray(y))

    def consensus(self, std_dev=3.0):
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
    x = asarray(x)
    y = asarray(y)

    xavg = x.mean(0)
    yavg = y.mean(0)

    err = y.std(0) * std_dev / sqrt(x.shape[0])

    return (xavg, yavg - err, yavg, yavg + err)


# def _consensus_curve(x, y, xlim=None, nparts=100, same_steps=False):
#     std_times = 3
#
#     if same_steps:
#         return _consensus_curve_same_steps(x, y, std_times)
#
#     x = np.asarray(x)
#     y = np.asarray(y)
#
#     if xlim:
#         xmin, xmax = xlim
#     else:
#         xmin, xmax = x.min(), x.max()
#
#     xstep = (xmax - xmin) / nparts
#
#     xavg = []
#     yavg = []
#
#     ytop = []
#     ybottom = []
#
#     for i in range(nparts):
#         left = xmin + i * xstep
#         right = left + xstep
#
#         xbulk, ybulk = _2dbulk(x, left, right, y)
#
#         if len(xbulk) > 0:
#             xavg.append(xbulk.mean())
#             yavg.append(ybulk.mean())
#
#             err = ybulk.std() * std_times / np.sqrt(len(ybulk))
#             ytop.append(yavg[-1] + err)
#             ybottom.append(yavg[-1] - err)
#
#     if xavg[0] > x.min():
#         xavg.insert(0, x.min())
#         yavg.insert(0, yavg[0])
#         ytop.insert(0, ytop[0])
#         ybottom.insert(0, ybottom[0])
#
#     if xavg[-1] < x.max():
#         xavg.append(x.max())
#         yavg.append(yavg[-1])
#         ytop.append(ytop[-1])
#         ybottom.append(ybottom[-1])
#
#     xavg = np.asarray(xavg)
#     yavg = np.asarray(yavg)
#
#     ytop = np.asarray(ytop)
#     ybottom = np.asarray(ybottom)
#
#     return (xavg, ybottom, yavg, ytop)
