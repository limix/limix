from numpy import asarray, concatenate, sort, zeros
from numpy.random import RandomState
from numpy.testing import assert_allclose
from pandas import DataFrame

from limix.stats import confusion_matrix


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
    (fpr, tpr) = cm.roc()
    assert_allclose(fpr[90], 0.729508196721)
    assert_allclose(tpr[90], 0.666666666667)
