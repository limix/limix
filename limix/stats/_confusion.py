from __future__ import division

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

    return ConfusionMatrix(P, N, ld_causal_markers, argsort(pv))


def getter(func):
    class ItemGetter(object):
        def __getitem__(self, i):
            return func(i)

        def __lt__(self, other):
            from numpy import s_

            return func(s_[:]) < other

        def __le__(self, other):
            from numpy import s_

            return func(s_[:]) <= other

        def __gt__(self, other):
            from numpy import s_

            return func(s_[:]) > other

        def __ge__(self, other):
            from numpy import s_

            return func(s_[:]) >= other

        def __eq__(self, other):
            from numpy import s_

            return func(s_[:]) == other

    return ItemGetter()


# TODO: document it
class ConfusionMatrix(object):
    def __init__(self, P, N, true_set, idx_rank):
        from numpy import empty, asarray, searchsorted

        self._TP = empty(P + N + 1, dtype=int)
        self._FP = empty(P + N + 1, dtype=int)
        if len(idx_rank) != P + N:
            raise ValueError(
                "Rank indices array has to have length equal" + " to ``P + N``."
            )

        true_set = asarray(true_set, int)
        true_set.sort()

        idx_rank = asarray(idx_rank, int)

        ins_pos = searchsorted(true_set, idx_rank)
        _confusion_matrix_tp_fp(P + N, ins_pos, true_set, idx_rank, self._TP, self._FP)
        self._N = N
        self._P = P

    @property
    def TP(self):
        return getter(lambda i: self._TP[i])

    @property
    def FP(self):
        return getter(lambda i: self._FP[i])

    @property
    def TN(self):
        return getter(lambda i: self._N - self.FP[i])

    @property
    def FN(self):
        return getter(lambda i: self._P - self.TP[i])

    @property
    def sensitivity(self):
        """ Sensitivity (also known as true positive rate.)
        """
        return getter(lambda i: self.TP[i] / self._P)

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
        return getter(lambda i: self.TN[i] / self._N)

    @property
    def precision(self):
        from numpy import nan

        return getter(
            lambda i: nan if i == 0 else self.TP[i] / (self.TP[i] + self.FP[i])
        )

    @property
    def npv(self):
        """ Negative predictive value.
        """
        return getter(lambda i: self.TN[i] / (self.TN[i] + self.FN[i]))

    @property
    def fallout(self):
        """ Fall-out (also known as false positive rate.)
        """
        return getter(lambda i: 1 - self.specifity[i])

    @property
    def fpr(self):
        return self.fallout

    @property
    def fnr(self):
        """ False negative rate.
        """
        return getter(lambda i: 1 - self.sensitivity[i])

    @property
    def fdr(self):
        """ False discovery rate.
        """
        return getter(lambda i: 1 - self.precision[i])

    @property
    def accuracy(self):
        """ Accuracy.
        """
        return getter(lambda i: (self.TP[i] + self.TN[i]) / (self._N + self._P))

    @property
    def f1score(self):
        """ F1 score (harmonic mean of precision and sensitivity).
        """

        def denominator(i):
            return 2 * self.TP[i] + self.FP[i] + self.FN[i]

        return getter(lambda i: 2 * self.TP[i] / denominator(i))

    def roc(self):
        from numpy import argsort

        tpr = self.tpr[1:]
        fpr = self.fpr[1:]

        idx = argsort(fpr)
        fpr = fpr[idx]
        tpr = tpr[idx]

        return (fpr, tpr)


def auc(fpr, tpr):
    left = fpr[0]
    area = 0.0
    for i in range(1, len(fpr)):
        width = fpr[i] - left
        area += width * tpr[i - 1]
        left = fpr[i]
    area += (1 - left) * tpr[-1]
    return area


def _confusion_matrix_tp_fp(n, ins_pos, true_set, idx_rank, TP, FP):
    TP[0] = 0
    FP[0] = 0
    i = 0
    while i < n:
        FP[i + 1] = FP[i]
        TP[i + 1] = TP[i]

        j = ins_pos[i]
        if j == len(true_set) or true_set[j] != idx_rank[i]:
            FP[i + 1] += 1
        else:
            TP[i + 1] += 1
        i += 1
