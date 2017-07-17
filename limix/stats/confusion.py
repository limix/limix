from __future__ import division

import logging

import numpy as np
from numba import jit
from numpy import argsort, asarray, concatenate, log10, nan, where
from numpy_sugar import is_crescent


@jit
def first_occurrence(arr, v):
    for i in range(arr.shape[0]):
        if arr[i] == v:
            return i
    return None


@jit
def _walk_left(pos, c, dist):
    assert dist > 0
    step = 0
    middle = pos[c]
    i = c
    while i > 0 and step < dist:
        i -= 1
        assert pos[i] <= middle
        step = (middle - pos[i])
    if step > dist:
        i += 1
    assert i <= c
    return i


@jit
def _walk_right(pos, c, dist):
    assert dist > 0
    step = 0
    middle = pos[c]
    i = c
    n = len(pos)
    while i < n - 1 and step < dist:
        i += 1
        assert pos[i] >= middle
        step = (pos[i] - middle)
    if step > dist:
        i -= 1
    assert i >= c
    return i


def roc_curve(multi_score, method, max_fpr=0.05):
    max_fpr = float(max_fpr)
    fprs = np.arange(0., max_fpr, step=0.001)
    tprs = np.empty_like(fprs)
    tprs_stde = np.empty_like(fprs)
    for (i, fpr) in enumerate(fprs):
        tprs_ = multi_score.get_tprs(method, fpr=fpr, approach='rank')
        tprs[i] = np.mean(tprs_)
        assert tprs[i] >= tprs[max(i - 1, 0)]
        tprs_stde[i] = np.std(tprs_) / np.sqrt(len(tprs_))
    return (fprs, tprs, tprs_stde)


def confusion_matrix(df, wsize=50000):
    """Provide a couple of scores based on the idea of windows around
       genetic markers.

       :param causals: Indices defining the causal markers.
       :param pos: Within-chromossome base-pair position of each candidate
                   marker, in crescent order.
    """
    logger = logging.getLogger(__name__)
    wsize = int(wsize)

    df.sort_values(by=['chrom', 'pos'], inplace=True)

    chromids = df['chrom'].unique()

    if len(chromids) == 0:
        raise ValueError("At least one chromossome is needed.")

    offset = 0
    idx_offset = 0
    pos = []
    causal = []
    pv = []
    for cid in sorted(chromids):
        df0 = df.query("chrom=='%s'" % cid)

        pos.append(offset + asarray(df0['pos'], float))
        pv.append(asarray(df0['pv'], float))
        offset += pos[-1][-1] + wsize // 2 + 1

        if df0['causal'].sum() > 0:
            causal.append(idx_offset + where(df0['causal'])[0])
            idx_offset += len(df0)

    pos = concatenate(pos)
    pv = concatenate(pv)
    causal = concatenate(causal)
    causal = asarray(causal, int)

    total_size = pos[-1] - pos[0]
    if wsize > 0.1 * total_size:
        perc = wsize // total_size * 100
        self._logger.warn(
            'The window size is %d%% of the total candidate' + ' region.',
            int(perc))

    ld_causal_markers = set()
    for (j, c) in enumerate(causal):
        if wsize == 1:
            right = left = pos[c]
        else:
            left = _walk_left(pos, c, wsize // 2)
            right = _walk_right(pos, c, wsize // 2)
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
            return func(np.s_[:]) < other

        def __le__(self, other):
            return func(np.s_[:]) <= other

        def __gt__(self, other):
            return func(np.s_[:]) > other

        def __ge__(self, other):
            return func(np.s_[:]) >= other

        def __eq__(self, other):
            return func(np.s_[:]) == other

    return ItemGetter()


class ConfusionMatrix(object):
    def __init__(self, P, N, true_set, idx_rank):
        self._TP = np.empty(P + N + 1, dtype=int)
        self._FP = np.empty(P + N + 1, dtype=int)
        assert len(idx_rank) == P + N

        true_set = np.asarray(true_set, int)
        true_set.sort()

        idx_rank = np.asarray(idx_rank, int)

        ins_pos = np.searchsorted(true_set, idx_rank)
        _confusion_matrix_tp_fp(P + N, ins_pos, true_set, idx_rank, self._TP,
                                self._FP)
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
        return getter(
            lambda i: (self.TP[i] + self.TN[i]) / (self._N + self._P))

    @property
    def f1score(self):
        """ F1 score (harmonic mean of precision and sensitivity).
        """
        return getter(
            lambda i: 2 * self.TP[i] / (2 * self.TP[i] + self.FP[i] + self.FN[i])
        )

    def roc(self):
        tpr = self.tpr[1:]
        fpr = self.fpr[1:]

        idx = np.argsort(fpr)
        fpr = fpr[idx]
        tpr = tpr[idx]

        return (fpr, tpr)


def auc(fpr, tpr):
    left = fpr[0]
    area = 0.
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


# if __name__ == '__main__':
#     import pandas as pd
#     import numpy as np
#
#     random = np.random.RandomState(0)
#
#     pos1 = np.sort(np.asarray(random.randint(1, 10000, size=100), float))
#     causal1 = np.zeros(100, dtype=bool)
#     causal1[1] = True
#     causal1[50] = True
#     pv1 = random.rand(100)
#
#     pos2 = np.sort(np.asarray(random.randint(1, 10000, size=25), float))
#     causal2 = np.zeros(25, dtype=bool)
#     causal2[3] = True
#     pv2 = random.rand(25)
#
#     df = pd.DataFrame(data=dict(
#         chrom=['1'] * 100 + ['2'] * 25,
#         pv=np.concatenate([pv1, pv2]),
#         pos=np.concatenate([pos1, pos2]),
#         causal=np.concatenate([causal1, causal2])))
#
#     cm = confusion_matrix(df, wsize=5)
#
#     print(cm.TP[90])
#     print(cm.FP[90])
#     print(cm.TN[90])
#     print(cm.FN[90])
#     print(cm.sensitivity[90])
#     print(cm.tpr[90])
#     print(cm.recall[90])
#     print(cm.specifity[90])
#     print(cm.precision[90])
#     print(cm.npv[90])
#     print(cm.fallout[90])
#     print(cm.fpr[90])
#     print(cm.fnr[90])
#     print(cm.fdr[90])
#     print(cm.accuracy[90])
#     print(cm.f1score[90])
#     (fpr, tpr) = cm.roc()
#     print(fpr[90])
#     print(tpr[90])
#     2
#     88
#     34
#     1
#     0.666666666667
#     0.666666666667
#     0.666666666667
#     0.27868852459
#     0.0222222222222
#     0.971428571429
#     0.72131147541
#     0.72131147541
#     0.333333333333
#     0.977777777778
#     0.288
#     0.0430107526882
#     0.729508196721
#     0.666666666667
