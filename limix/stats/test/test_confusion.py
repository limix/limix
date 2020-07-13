import pytest
from numpy import asarray, concatenate, sort, zeros, nan, inf
from numpy.random import RandomState
from numpy.testing import assert_allclose
from pandas import DataFrame

from limix.stats import confusion_matrix, ConfusionMatrix


def test_stats_confusion():
    random = RandomState(0)

    pos1 = sort(asarray(random.randint(1, 10000, size=100), float))
    causal1 = zeros(100, dtype=bool)
    causal1[1] = True
    causal1[50] = True
    pv1 = random.rand(100)

    pos2 = sort(asarray(random.randint(1, 10000, size=25), float))
    causal2 = zeros(25, dtype=bool)
    causal2[3] = True
    pv2 = random.rand(25)

    df = DataFrame(
        data=dict(
            chrom=["1"] * 100 + ["2"] * 25,
            pv=concatenate([pv1, pv2]),
            pos=concatenate([pos1, pos2]),
            causal=concatenate([causal1, causal2]),
        )
    )

    cm = confusion_matrix(df, wsize=5)

    assert_allclose(cm.precision[0], nan)
    assert_allclose(cm.npv[-1], nan)

    assert_allclose(cm.TP[0], 0)
    assert_allclose(cm.FP[0], 0)
    assert_allclose(cm.FN[0], 3)
    assert_allclose(cm.TN[0], 122)

    assert_allclose(cm.TP[-1], 3)
    assert_allclose(cm.FP[-1], 122)
    assert_allclose(cm.FN[-1], 0)
    assert_allclose(cm.TN[-1], 0)

    assert_allclose(cm.TP[90], 2)
    assert_allclose(cm.FP[90], 88)
    assert_allclose(cm.TN[90], 34)
    assert_allclose(cm.FN[90], 1)
    assert_allclose(cm.sensitivity[90], 0.666666666667)
    assert_allclose(cm.tpr[90], 0.666666666667)
    assert_allclose(cm.recall[90], 0.666666666667)
    assert_allclose(cm.specifity[90], 0.27868852459)
    assert_allclose(cm.precision[90], 0.0222222222222)
    assert_allclose(cm.npv[90], 0.971428571429)
    assert_allclose(cm.fallout[90], 0.72131147541)
    assert_allclose(cm.fpr[90], 0.72131147541)
    assert_allclose(cm.fnr[90], 0.333333333333)
    assert_allclose(cm.fdr[90], 0.977777777778)
    assert_allclose(cm.accuracy[90], 0.288)
    assert_allclose(cm.f1score[90], 0.0430107526882)
    roc = cm.roc()
    assert_allclose(roc.fpr[90], 0.729508196721)
    assert_allclose(roc.tpr[90], 0.666666666667)

    assert_allclose(roc.auc, 0.28415300546448086)


def test_stats_confusion_matrix():
    true_samples = [3, 2]
    N = 100
    sorted_samples = [1, 3, 100]

    cm = ConfusionMatrix(true_samples, N, sorted_samples)

    assert_allclose(cm.precision[0], nan)
    assert_allclose(cm.precision[-1], 1 / 3)

    assert_allclose(cm.npv[0], 0.9803921568627451)
    assert_allclose(cm.npv[-1], nan)

    assert_allclose(cm.sensitivity[0], 0.0)
    assert_allclose(cm.sensitivity[-1], 0.5)

    assert_allclose(cm.specifity[0], 1.0)
    assert_allclose(cm.specifity[-1], 0.98)

    assert_allclose(cm.f1score[0], 0.0)
    assert_allclose(cm.f1score[-1], 0.4)

    assert_allclose(cm.accuracy[0], 0.9803921568627451)
    assert_allclose(cm.accuracy[-1], 0.9705882352941176)

    roc = cm.roc()
    assert_allclose(roc.fpr, [0.01000000, 0.01000000, 0.02000000])
    assert_allclose(roc.tpr, [0.00000000, 0.50000000, 0.50000000])

    assert_allclose(roc.auc, 0.495)


def test_stats_confusion_matrix_len1():
    true_samples = [3, 2]
    N = 100
    sorted_samples = [1]

    cm = ConfusionMatrix(true_samples, N, sorted_samples)

    assert_allclose(cm.precision[0], nan)
    assert_allclose(cm.precision[-1], 0.0)

    assert_allclose(cm.npv[0], 0.9803921568627451)
    assert_allclose(cm.npv[-1], nan)

    assert_allclose(cm.sensitivity[0], 0.0)
    assert_allclose(cm.sensitivity[-1], 0.0)

    assert_allclose(cm.specifity[0], 1.0)
    assert_allclose(cm.specifity[-1], 0.99)

    assert_allclose(cm.f1score[0], 0.0)
    assert_allclose(cm.f1score[-1], 0.0)

    assert_allclose(cm.accuracy[0], 0.9803921568627451)
    assert_allclose(cm.accuracy[-1], 0.9705882352941176)

    roc = cm.roc()
    assert_allclose(roc.fpr, [0.01000000])
    assert_allclose(roc.tpr, [0.00000000])

    assert_allclose(roc.auc, 0.0)


def test_stats_confusion_matrix_limit():
    true_samples = [3, 2]
    sorted_samples = [1]

    with pytest.raises(ValueError):
        ConfusionMatrix(true_samples, 0, sorted_samples)

    cm = ConfusionMatrix(true_samples, 0, [])

    assert_allclose(cm.precision[0], nan)

    assert_allclose(cm.npv[0], nan)

    assert_allclose(cm.sensitivity[0], 0.0)

    assert_allclose(cm.specifity[0], nan)

    assert_allclose(cm.f1score[0], 0.0)

    assert_allclose(cm.accuracy[0], 0.0)

    with pytest.raises(ValueError):
        cm.roc()
