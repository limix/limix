from limix.glmm import GLMMComposer
from limix.qc import normalise_covariance
from numpy import dot, ones, eye
from numpy.testing import assert_allclose
from numpy.random import RandomState


def test_glmm_composer():
    random = RandomState(0)
    nsamples = 50

    glmm = GLMMComposer(nsamples)

    glmm.fixed_effects.append_offset()

    X0 = random.randn(nsamples)
    glmm.fixed_effects.append(X0)
    glmm.fixed_effects[0].offset = 1
    glmm.fixed_effects[1].effsizes = [1]

    assert_allclose(glmm.fixed_effects.mean.value() - X0, ones(nsamples))

    X12 = random.randn(nsamples, 2)
    glmm.fixed_effects.append(X12)

    G0 = random.randn(nsamples, 100)
    K0 = normalise_covariance(dot(G0, G0.T))
    glmm.covariance_matrices.append(K0)

    G1 = random.randn(nsamples, 100)
    K1 = normalise_covariance(dot(G1, G1.T))
    glmm.covariance_matrices.append(K1)

    glmm.covariance_matrices.append_iid_noise()
    glmm.covariance_matrices[0].scale = 1
    glmm.covariance_matrices[1].scale = 0
    glmm.covariance_matrices[2].scale = 1
    K = glmm.covariance_matrices.cov.value()
    assert_allclose(K, K0 + eye(nsamples))

    y = random.randn(nsamples)
    glmm.y = y

    glmm.fit(verbose=False)

    assert_allclose(glmm.covariance_matrices[0].scale, 0, atol=1e-6)
    assert_allclose(glmm.covariance_matrices[1].scale, 0, atol=1e-6)
    assert_allclose(glmm.covariance_matrices[2].scale, 1.099905167170892, atol=1e-6)

    assert_allclose(glmm.lml(), -73.32753446649403, atol=1e-6)


# def test_glmm_composer_plot():
#     random = RandomState(0)
#     nsamples = 50

#     glmm = GLMMComposer(nsamples)

#     glmm.fixed_effects.append_offset()

#     X0 = random.randn(nsamples)
#     glmm.fixed_effects.append(X0)
#     glmm.fixed_effects[0].offset = 1
#     glmm.fixed_effects[1].effsizes = [1]

#     assert_allclose(glmm.fixed_effects.mean.value() - X0, ones(nsamples))

#     X12 = random.randn(nsamples, 2)
#     glmm.fixed_effects.append(X12)

#     G0 = random.randn(nsamples, 100)
#     K0 = normalise_covariance(dot(G0, G0.T))
#     glmm.covariance_matrices.append(K0)

#     G1 = random.randn(nsamples, 100)
#     K1 = normalise_covariance(dot(G1, G1.T))
#     glmm.covariance_matrices.append(K1)

#     glmm.covariance_matrices.append_iid_noise()
#     glmm.covariance_matrices[0].scale = 1
#     glmm.covariance_matrices[1].scale = 0
#     glmm.covariance_matrices[2].scale = 1
#     K = glmm.covariance_matrices.cov.value()
#     assert_allclose(K, K0 + eye(nsamples))

#     y = random.randn(nsamples)
#     glmm.y = y

#     glmm.fit(verbose=True)

#     assert_allclose(glmm.covariance_matrices[0].scale, 0, atol=1e-6)
#     assert_allclose(glmm.covariance_matrices[1].scale, 0, atol=1e-6)
#     assert_allclose(glmm.covariance_matrices[2].scale, 1.099905167170892, atol=1e-6)

#     assert_allclose(glmm.lml(), -73.32753446649403, atol=1e-6)

#     print()
#     print(glmm)
#     # glmm.plot()
#     # print(glmm.decomp())
