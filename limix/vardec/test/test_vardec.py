from numpy import concatenate, exp, eye, ones, zeros
from numpy.random import RandomState
from numpy.testing import assert_allclose

from limix.stats import multivariate_normal as mvn
from limix.vardec import VarDec


def test_vardec():
    random = RandomState(0)
    nsamples = 20

    X = random.randn(nsamples, 2)
    X = (X - X.mean(0)) / X.std(0)
    X = concatenate((ones((nsamples, 1)), X), axis=1)
    lik = "normal"

    K0 = random.randn(nsamples, 10)
    K0 = K0 @ K0.T
    K0 /= K0.diagonal().mean()
    K0 += eye(nsamples) * 1e-4

    K1 = random.randn(nsamples, 10)
    K1 = K1 @ K1.T
    K1 /= K1.diagonal().mean()
    K1 += eye(nsamples) * 1e-4

    y = X @ random.randn(3) + mvn(random, zeros(nsamples), K0)
    y += mvn(random, zeros(nsamples), K1)

    vardec = VarDec(y, lik, X)
    vardec.append(K0)
    vardec.append(K1)
    vardec.append_iid()

    vardec.fit(verbose=False)
    assert_allclose(vardec.covariance[0].scale, 0.38473522809891697)
    assert_allclose(vardec.covariance[1].scale, 1.1839796169221422)
    assert_allclose(vardec.covariance[2].scale, 2.061153622438558e-09, atol=1e-5)
    assert_allclose(vardec.lml(), -21.91827344966165)
    assert_allclose(
        vardec.effsizes, [-0.5008873352111712, -1.193536406235688, -0.28254298530554534]
    )


def test_vardec_2_matrices():
    random = RandomState(0)
    nsamples = 20

    X = random.randn(nsamples, 2)
    X = (X - X.mean(0)) / X.std(0)
    X = concatenate((ones((nsamples, 1)), X), axis=1)
    lik = "normal"

    K = random.randn(nsamples, 10)
    K = K @ K.T
    K /= K.diagonal().mean()
    K += eye(nsamples) * 1e-4

    y = X @ random.randn(3) + mvn(random, zeros(nsamples), K) + random.randn(nsamples)

    vardec = VarDec(y, lik, X)
    vardec.append(K)
    vardec.append_iid()

    vardec.fit(verbose=False)
    assert_allclose(vardec.covariance[0].scale, 0.32199728815536727, rtol=1e-5)
    assert_allclose(vardec.covariance[1].scale, 1.4182987383374532, rtol=1e-5)
    assert_allclose(vardec.lml(), -33.63946372828994, rtol=1e-5)


def test_vardec_poisson():
    random = RandomState(0)
    nsamples = 20

    X = random.randn(nsamples, 2)
    X = (X - X.mean(0)) / X.std(0)
    X = concatenate((ones((nsamples, 1)), X), axis=1)
    lik = "poisson"

    K0 = random.randn(nsamples, 10)
    K0 = K0 @ K0.T
    K0 /= K0.diagonal().mean()
    K0 += eye(nsamples) * 1e-4

    K1 = random.randn(nsamples, 10)
    K1 = K1 @ K1.T
    K1 /= K1.diagonal().mean()
    K1 += eye(nsamples) * 1e-4

    y = X @ random.randn(3) + mvn(random, zeros(nsamples), K0)
    y += mvn(random, zeros(nsamples), K1)
    y = exp((y - y.mean()) / y.std())

    vardec = VarDec(y, lik, X)
    vardec.append(K0)
    vardec.append(K1)
    vardec.append_iid()

    vardec.fit(verbose=False)
    assert_allclose(vardec.covariance[0].scale, 2.808478303826397e-09, atol=1e-5)
    assert_allclose(vardec.covariance[1].scale, 0.3503800007985209)
    assert_allclose(vardec.covariance[2].scale, 2.061153622438558e-09, atol=1e-5)
    assert_allclose(vardec.lml(), -28.887285796984564)


def test_vardec_poisson_2_matrices():
    random = RandomState(0)
    nsamples = 20

    X = random.randn(nsamples, 2)
    X = (X - X.mean(0)) / X.std(0)
    X = concatenate((ones((nsamples, 1)), X), axis=1)
    lik = "poisson"

    K = random.randn(nsamples, 10)
    K = K @ K.T
    K /= K.diagonal().mean()
    K += eye(nsamples) * 1e-4

    y = X @ random.randn(3) + mvn(random, zeros(nsamples), K)
    y = exp((y - y.mean()) / y.std())

    vardec = VarDec(y, lik, X)
    vardec.append(K)
    vardec.append_iid()

    vardec.fit(verbose=False)
    assert_allclose(vardec.covariance[0].scale, 0.10726852397002325, atol=1e-5)
    assert_allclose(vardec.covariance[1].scale, 4.168569936272955e-11, atol=1e-5)
    assert_allclose(vardec.lml(), -26.36419072811823)
