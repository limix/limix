from limix.vardec import VarDec
from limix.qc import normalise_covariance
from numpy import dot, ones, eye, concatenate, zeros
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


# def test_glmm_composer_plot():
#     random = RandomState(0)
#     nsamples = 50

#     vardec = VarDec(nsamples)

#     vardec.fixed_effects.append_offset()

#     X0 = random.randn(nsamples)
#     vardec.fixed_effects.append(X0)
#     vardec.fixed_effects[0].offset = 1
#     vardec.fixed_effects[1].effsizes = [1]

#     assert_allclose(vardec.fixed_effects.mean.value() - X0, ones(nsamples))

#     X12 = random.randn(nsamples, 2)
#     vardec.fixed_effects.append(X12)

#     G0 = random.randn(nsamples, 100)
#     K0 = normalise_covariance(dot(G0, G0.T))
#     vardec.covariance_matrices.append(K0)

#     G1 = random.randn(nsamples, 100)
#     K1 = normalise_covariance(dot(G1, G1.T))
#     vardec.covariance_matrices.append(K1)

#     vardec.covariance_matrices.append_iid_noise()
#     vardec.covariance_matrices[0].scale = 1
#     vardec.covariance_matrices[1].scale = 0
#     vardec.covariance_matrices[2].scale = 1
#     K = vardec.covariance_matrices.cov.value()
#     assert_allclose(K, K0 + eye(nsamples))

#     y = random.randn(nsamples)
#     vardec.y = y

#     vardec.fit(verbose=True)

#     assert_allclose(vardec.covariance_matrices[0].scale, 0, atol=1e-6)
#     assert_allclose(vardec.covariance_matrices[1].scale, 0, atol=1e-6)
#     assert_allclose(vardec.covariance_matrices[2].scale, 1.099905167170892, atol=1e-6)

#     assert_allclose(vardec.lml(), -73.32753446649403, atol=1e-6)

#     print()
#     print(vardec)
#     # vardec.plot()
#     # print(vardec.decomp())
