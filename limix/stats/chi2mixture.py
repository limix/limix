import pdb

import scipy as SP
import scipy.stats as STATS


class Chi2mixture(object):
    def __init__(self,
                 scale_min=0.1,
                 scale_max=5.0,
                 dof_min=0.1,
                 dof_max=5.0,
                 n_intervals=100,
                 qmax=0.1,
                 tol=0):
        """
        class for a mixture of chi2-distributed random variables of the form

        (1-mixture)*chi2(0) + (mixture)*scale*chi2(dof),

        where
        scale   is the scaling parameter for the scales chi-square distribution
        dof     are the degrees of freedom of the second component
        mixture is the probability of beeing in the first component

        input:
        scale_min   minimum value used for fitting the scale parameter
        scale_max   maximum value used for fitting the scale parameter
        dofmin      minimum value used for fitting the dof parameter
        dofmax      maximum value used for fitting the dof parameter
        qmax        only the top qmax quantile is used for the fit
        n_interval  number of intervals when performing gridsearch
        tol        tolerance of being zero
        """

        self.scale_min = scale_min
        self.scale_max = scale_max
        self.dof_min = dof_min
        self.dof_max = dof_max
        self.qmax = qmax
        self.n_intervals = n_intervals
        self.tol = tol

    def estimate_chi2mixture(self, lrt):
        """
        estimates the parameters of a mixture of a chi-squared random variable of degree
        0 and a scaled chi-squared random variable of degree d

        (1-mixture)*chi2(0) + (mixture)*scale*chi2(dof),

        where
        scale   is the scaling parameter for the scales chi-square distribution
        dof     are the degrees of freedom of the second component
        mixture is the probability of beeing in the first component

        input:
        lrt         [Ntests] vector of test statistics
        """
        """
        step 1: estimate the probability of being in component one
        """
        self.mixture = 1 - (lrt <= self.tol).mean()
        n_false = SP.sum(lrt > self.tol)
        """
        step 2: only use the largest qmax fraction of test statistics to estimate the
                remaining parameters
        """
        n_fitting = SP.ceil(self.qmax * n_false)
        lrt_sorted = -SP.sort(-lrt)[:n_fitting]
        q = SP.linspace(0, 1, n_false)[1:n_fitting + 1]
        log_q = SP.log10(q)
        """
        step 3: fitting scale and dof by minimizing the squared error of the log10 p-values
                with their theorietical values [uniform distribution]
        """
        MSE_opt = SP.inf
        MSE = SP.zeros((self.n_intervals, self.n_intervals))

        for i, scale in enumerate(
                SP.linspace(self.scale_min, self.scale_max, self.n_intervals)):
            for j, dof in enumerate(
                    SP.linspace(self.dof_min, self.dof_max, self.n_intervals)):
                p = STATS.chi2.sf(lrt_sorted / scale, dof)
                log_p = SP.log10(p)
                MSE[i, j] = SP.mean((log_q - log_p)**2)
                if MSE[i, j] < MSE_opt:
                    MSE_opt = MSE[i, j]
                    self.scale = scale
                    self.dof = dof

    def sf(self, lrt):
        """
        computes the survival function of a mixture of a chi-squared random variable of degree
        0 and a scaled chi-squared random variable of degree d
        """
        _lrt = SP.copy(lrt)
        _lrt[lrt < self.tol] = 0
        pv = self.mixture * STATS.chi2.sf(_lrt / self.scale, self.dof)
        return pv


if __name__ == "__main__":
    scale = 0.3
    dof = 2
    mixture = 0.2
    n_test = 100000
    n_train = 100000

    lrt_train = SP.zeros((n_train))
    lrt_test = SP.zeros((n_test))

    for i in range(dof):
        x = SP.random.randn(n_train)
        lrt_train += scale * (x * x)
        x = SP.random.randn(n_test)
        lrt_test += scale * (x * x)

    idx_chi2_0 = SP.random.permutation(n_train)[:(1 - mixture) * n_train]
    lrt_train[idx_chi2_0] = 1e-10 * SP.random.randn((1 - mixture) * n_train)

    idx_chi2_0 = SP.random.permutation(n_test)[:(1 - mixture) * n_test]
    lrt_test[idx_chi2_0] = 1e-10 * SP.random.randn((1 - mixture) * n_test)

    chi2mix = Chi2mixture(
        scale_min=0.1,
        scale_max=5.0,
        dof_min=0.1,
        dof_max=5.0,
        qmax=0.1,
        tol=4e-3)
    chi2mix.estimate_chi2mixture(lrt_train)
    chi2mix.sf(lrt_test)

    print(('true scale = %.2f' % scale))
    print(('true dof   = %.2f' % dof))
    print(('true mixt  = %.2f' % mixture))
    print(('est scale  = %.2f' % chi2mix.scale))
    print(('est dof    = %.2f' % chi2mix.dof))
    print(('est mixt   = %.2f' % chi2mix.mixture))

    import matplotlib.pylab as PLT
    fig = PLT.figure()
    fig.add_subplot(121)
    PLT.hist(chi2mix.sf(lrt_train), normed=True)
    fig.add_subplot(122)
    PLT.hist(chi2mix.sf(lrt_test), normed=True)

    pval_test = chi2mix.sf(lrt_test)
    print('print Type 1 error estimate')
    alpha = 1e-5
    print('... alpha=1e-5')
    print(('...... type 1 error: %.2e' % (pval_test < alpha).mean()))
    alpha = 1e-4
    print('... alpha=1e-4')
    print(('...... type 1 error: %.2e' % (pval_test < alpha).mean()))
    alpha = 1e-3
    print('... alpha=1e-3')
    print(('...... type 1 error: %.2e' % (pval_test < alpha).mean()))
