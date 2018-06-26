from limix.glmm import GLMMComposer
from numpy import dot, sqrt, zeros
from numpy.random import RandomState


def test_glmm_composer():
    random = RandomState(0)
    nsamples = 50

    glmm = GLMMComposer()

    glmm.fixed_effects.append_offset()

    X0 = random.randn(nsamples)
    glmm.fixed_effects.append(X0)

    X12 = random.randn(nsamples, 2)
    glmm.fixed_effects.append(X12)

    G0 = random.randn(nsamples, 100)
    K0 = dot(G0, G0.T)
    glmm.covariance_matrices.append(K0)

    G1 = random.randn(nsamples, 100)
    K1 = dot(G1, G1.T)
    glmm.covariance_matrices.append(K1)

    glmm.covariance_matrices.append_iid_noise()

    y = random.randn(nsamples)
    glmm.y = y
