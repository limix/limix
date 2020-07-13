from __future__ import division
from __future__ import annotations

from typing import Iterable

import logging


def _get_jit():
    try:
        from numba import jit
    except ImportError:

        def jit(x, *args, **kwargs):
            return x

    return jit


def _get_walk_left():
    jit = _get_jit()

    @jit(cache=True)
    def _walk_left(pos, c, dist):
        step = 0
        middle = pos[c]
        i = c
        while i > 0 and step < dist:
            i -= 1
            step = middle - pos[i]
        if step > dist:
            i += 1
        return i

    return _walk_left


def _get_walk_right():
    jit = _get_jit()

    @jit(cache=True)
    def _walk_right(pos, c, dist):
        step = 0
        middle = pos[c]
        i = c
        n = len(pos)
        while i < n - 1 and step < dist:
            i += 1
            step = pos[i] - middle
        if step > dist:
            i -= 1
        return i

    return _walk_right


def roc_curve(multi_score, method, max_fpr=0.05):
    from numpy import arange, empty_like, mean, std, sqrt

    max_fpr = float(max_fpr)
    fprs = arange(0.0, max_fpr, step=0.001)
    tprs = empty_like(fprs)
    tprs_stde = empty_like(fprs)
    for (i, fpr) in enumerate(fprs):
        tprs_ = multi_score.get_tprs(method, fpr=fpr, approach="rank")
        tprs[i] = mean(tprs_)
        tprs_stde[i] = std(tprs_) / sqrt(len(tprs_))
    return (fprs, tprs, tprs_stde)


# TODO: convert to numpy style
def confusion_matrix(df, wsize=50000):
    """Provide a couple of scores based on the idea of windows around
       genetic markers.

       :param causals: Indices defining the causal markers.
       :param pos: Within-chromossome base-pair position of each candidate
                   marker, in crescent order.
    """
    from numpy import argsort, asarray, concatenate, where

    logger = logging.getLogger(__name__)
    wsize = int(wsize)

    if "chrom" not in df:
        df = df.assign(chrom=["1"] * len(df))

    df.sort_values(by=["chrom", "pos"], inplace=True)

    chromids = df["chrom"].unique()

    offset = 0
    idx_offset = 0
    pos = []
    causal = []
    pv = []
    for cid in sorted(chromids):
        df0 = df.query("chrom=='%s'" % cid)

        pos.append(offset + asarray(df0["pos"], float))
        pv.append(asarray(df0["pv"], float))
        offset += pos[-1][-1] + wsize // 2 + 1

        if df0["causal"].sum() > 0:
            causal.append(idx_offset + where(df0["causal"])[0])
            idx_offset += len(df0)

    pos = concatenate(pos)
    pv = concatenate(pv)
    causal = concatenate(causal)
    causal = asarray(causal, int)

    total_size = pos[-1] - pos[0]
    if wsize > 0.1 * total_size:
        perc = wsize // total_size * 100
        logger.warn(
            "The window size is %d%% of the total candidate" + " region.", int(perc)
        )

    ld_causal_markers = set()
    for _, c in enumerate(causal):
        if wsize == 1:
            right = left = pos[c]
        else:
            left = _get_walk_left()(pos, c, wsize // 2)
            right = _get_walk_right()(pos, c, wsize // 2)
        for i in range(left, right + 1):
            ld_causal_markers.add(i)

    P = len(ld_causal_markers)
    N = len(pos) - P

    ld_causal_markers = list(ld_causal_markers)

    logger.info("Found %d positive and %d negative markers.", P, N)

    return ConfusionMatrix(ld_causal_markers, N, argsort(pv))


class ROC:
    def __init__(self, fpr: Iterable[float], tpr: Iterable[float]):
        from numpy import asarray

        self._fpr = asarray(fpr, float)
        self._tpr = asarray(tpr, float)

    @property
    def fpr(self):
        return self._fpr

    @property
    def tpr(self):
        return self._tpr

    @property
    def auc(self) -> float:
        left = self.fpr[0]
        area = 0.0
        for i in range(1, len(self.fpr)):
            width = self.fpr[i] - left
            area += width * self.tpr[i - 1]
            left = self.fpr[i]
        area += (1 - left) * self.tpr[-1]
        return area


class ConfusionMatrix:
    """
    Confusion matrix.

    Parameters
    ----------
    true_samples
        Set of all positive samples from the solution space.
    N
        Number of negative samples.
    sorted_samples
        Samples sorted from the most to the least likely one to be considered positive.
    """

    def __init__(
        self, true_samples: Iterable[int], N: int, sorted_samples: Iterable[int]
    ):
        from numpy import empty, asarray

        if len(set(sorted_samples) - set(true_samples)) > N:
            raise ValueError("Invalid number of negative samples.")

        true_arr = asarray(true_samples, int)
        P = len(true_arr)

        sorted_arr = asarray(sorted_samples, int)

        self._TP = empty(len(sorted_arr) + 1, int)
        self._FP = empty(len(sorted_arr) + 1, int)

        self._N = N
        self._P = P
        self._num_sorted_samples = len(sorted_arr)
        self._set_tp_fp(true_arr, sorted_arr)

    def _set_tp_fp(self, true_samples, sorted_samples):
        from numpy import searchsorted, asarray

        true_arr = asarray(true_samples, int)
        true_arr.sort()

        sorted_arr = asarray(sorted_samples, int)
        ins_pos = searchsorted(true_arr, sorted_arr)

        self._TP[0] = 0
        self._FP[0] = 0
        i = 0
        while i < len(sorted_arr):
            self._FP[i + 1] = self._FP[i]
            self._TP[i + 1] = self._TP[i]

            j = ins_pos[i]
            if j == len(true_arr) or true_arr[j] != sorted_arr[i]:
                self._FP[i + 1] += 1
            else:
                self._TP[i + 1] += 1
            i += 1

    @property
    def TP(self):
        return self._TP

    @property
    def FP(self):
        return self._FP

    @property
    def TN(self):
        return self._N - self.FP

    @property
    def FN(self):
        return self._P - self.TP

    @property
    def sensitivity(self):
        """ Sensitivity (also known as true positive rate.)
        """
        return self.TP / self._P

    @property
    def tpr(self):
        return self.sensitivity

    @property
    def recall(self):
        return self.sensitivity

    @property
    def specifity(self):
        """ Specifity (also known as true negative rate.)
        """
        return self.TN / self._N

    @property
    def precision(self):
        from numpy import nan, empty

        r = empty(self._num_sorted_samples + 1)
        r[0] = nan
        r[1:] = self.TP[1:] / (self.TP[1:] + self.FP[1:])

        return r

    @property
    def npv(self):
        """ Negative predictive value.
        """
        from numpy import nan, empty

        r = empty(self._num_sorted_samples + 1)
        r[-1] = nan
        r[:-1] = self.TN[:-1] / (self.TN[:-1] + self.FN[:-1])

        return r

    @property
    def fallout(self):
        """ Fall-out (also known as false positive rate.)
        """
        return 1 - self.specifity

    @property
    def fpr(self):
        return self.fallout

    @property
    def fnr(self):
        """ False negative rate.
        """
        return 1 - self.sensitivity

    @property
    def fdr(self):
        """ False discovery rate.
        """
        return 1 - self.precision

    @property
    def accuracy(self):
        """ Accuracy.
        """
        return (self.TP + self.TN) / (self._N + self._P)

    @property
    def f1score(self):
        """ F1 score (harmonic mean of precision and sensitivity).
        """
        return 2 * self.TP / (2 * self.TP + self.FP + self.FN)

    def roc(self) -> ROC:
        from numpy import argsort

        if self._num_sorted_samples < 1:
            raise ValueError("Not enough sorted samples.")

        tpr = self.tpr[1:]
        fpr = self.fpr[1:]

        idx = argsort(fpr)
        fpr = fpr[idx]
        tpr = tpr[idx]

        return ROC(fpr, tpr)
