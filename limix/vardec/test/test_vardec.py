from limix.vardec import VarDec
from numpy import ones, eye, concatenate, zeros, exp
from numpy.testing import assert_allclose
from numpy.random import RandomState


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

    mvn = random.multivariate_normal
    y = X @ random.randn(3) + mvn(zeros(nsamples), K0) + mvn(zeros(nsamples), K1)

    vardec = VarDec(y, lik, X)
    vardec.append(K0)
    vardec.append(K1)
    vardec.append_iid()

    vardec.fit(verbose=False)
    assert_allclose(vardec.covariance[0].scale, 0.42493502300821745)
    assert_allclose(vardec.covariance[1].scale, 1.775872164537344)
    assert_allclose(vardec.covariance[2].scale, 2.061153622438558e-09, atol=1e-5)
    assert_allclose(vardec.lml(), -24.447408443017064)


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

    mvn = random.multivariate_normal
    y = X @ random.randn(3) + mvn(zeros(nsamples), K) + random.randn(nsamples)

    vardec = VarDec(y, lik, X)
    vardec.append(K)
    vardec.append_iid()

    vardec.fit(verbose=False)
    assert_allclose(vardec.covariance[0].scale, 0.5331692582164862, rtol=1e-5)
    assert_allclose(vardec.covariance[1].scale, 1.45673841962057, rtol=1e-5)
    assert_allclose(vardec.lml(), -34.694171078044846, rtol=1e-5)


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

    mvn = random.multivariate_normal
    y = X @ random.randn(3) + mvn(zeros(nsamples), K0) + mvn(zeros(nsamples), K1)
    y = exp((y - y.mean()) / y.std())

    vardec = VarDec(y, lik, X)
    vardec.append(K0)
    vardec.append(K1)
    vardec.append_iid()

    vardec.fit(verbose=False)
    assert_allclose(vardec.covariance[0].scale, 2.6905366173575983e-09, atol=1e-5)
    assert_allclose(vardec.covariance[1].scale, 0.3965071579047076)
    assert_allclose(vardec.covariance[2].scale, 2.061153622438558e-09, atol=1e-5)
    assert_allclose(vardec.lml(), -29.55133669213795)


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

    mvn = random.multivariate_normal
    y = X @ random.randn(3) + mvn(zeros(nsamples), K)
    y = exp((y - y.mean()) / y.std())

    vardec = VarDec(y, lik, X)
    vardec.append(K)
    vardec.append_iid()

    vardec.fit(verbose=False)
    assert_allclose(vardec.covariance[0].scale, 0.0009999999008707368, atol=1e-5)
    assert_allclose(vardec.covariance[1].scale, 9.912926347978695e-11, atol=1e-5)
    assert_allclose(vardec.lml(), -28.059529339760136)
