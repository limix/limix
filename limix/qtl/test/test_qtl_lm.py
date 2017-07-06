# from numpy import dot, sqrt, zeros
# from numpy.random import RandomState
# from numpy.testing import assert_allclose
#
# from limix.qtl import qtl_test_lm
#
#
# def test_qtl_lm():
#     random = RandomState(0)
#     nsamples = 50
#
#     G = random.randn(50, 100)
#
#     y = dot(G, random.randn(100)) / sqrt(100) + 0.2 * random.randn(nsamples)
#
#     M = G[:, :5]
#     X = G[:, 5:]
#     import pdb
#     pdb.set_trace()
#     lmm = qtl_test_lm(X, y, M=M, verbose=True)
#     # pv = lmm.variant_pvalues
#     # assert_allclose(pv, [0.41890603984302488, 0.71500332916140363], rtol=1e-3)
