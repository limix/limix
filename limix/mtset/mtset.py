import h5py
import scipy as sp
import scipy.linalg as la
import time as TIME
import warnings
import os


class MTSet():
    r"""
    Multi-trait set tests (mtSet).

    Mixed-model approach that enables joint analysis across multiple correlated
    traits while accounting for population structure and relatedness [CRS15]_.

    Parameters
    ----------
    Y : (`N`, `P`) ndarray
        Individuals by phenotypes.
    F : (`N`, `K`) ndarray, optional
        Individuals by covariates. By default, no fixed effect term is set.
    R : (`N`, `N`) ndarray, optional
        LMM-covariance/kinship coefficients.
        ``U_R`` and ``S_R`` can be provided instead of ``R``.
        If neither ``R`` nor ``U_R`` and ``S_R`` are provided,
        the null models has iid normal residuals.
    U_R : (`N`, `N`) ndarray, optional
        Eigenvectors of ``R``.
        ``U_R`` and ``S_R`` can be provided instead of ``R``.
        If neither ``R`` nor ``U_R`` and ``S_R`` are provided,
        iid normal residuals are considered.
    S_R : (`N`, ) ndarray, optional
        Eigenvalues of ``R``.
        ``U_R`` and ``S_R`` can be provided instead of ``R``.
        If neither ``R`` nor ``U_R`` and ``S_R`` are provided,
        iid normal residuals are considered.
    traitID : ndarray, optional
        Phenotype IDs.
    rank : int, optional
        Rank of the trait set covariance. The default value is 1.

    Examples
    --------
    This example shows how to fit mtSet when modelling population
    structure/relatedness using the full genetic relatedness matrix.

    .. doctest::

        >>> from numpy.random import RandomState
        >>> from limix.mtSet import MTSet
        >>> from numpy import dot, eye, ones, concatenate
        >>> from numpy import set_printoptions
        >>> set_printoptions(4)
        >>>
        >>> random = RandomState(1)
        >>>
        >>> N = 100
        >>> P = 2
        >>>
        >>> pheno = random.randn(N, P)
        >>> pheno-= pheno.mean(0)
        >>> pheno/= pheno.std(0)
        >>>
        >>> W = random.randn(N, 10)
        >>> kinship = dot(W, W.T) / float(10)
        >>> kinship+= 1e-4 * eye(N)
        >>>
        >>> mtset = MTSet(pheno, R=kinship)
        >>> res_null = mtset.fitNull()
        >>>
        >>> print(res_null['conv'][0])
        True
        >>> print('%.2f'%res_null['NLL0'])
        98.43
        >>>
        >>> S = 4
        >>> snp_set = (random.rand(N, S) < 0.2).astype(float)
        >>> res = mtset.optimize(snp_set)
        >>>
        >>> print("%.4f" % res['LLR'][0])
        0.0616
        >>> print(res['Cr'])
        [[ 0.0006  0.0025]
         [ 0.0025  0.0113]]

    This example shows how to fit mtSet when modelling population
    structure/relatedness by introducing the top principle components
    of the genetic relatedness matrix (``pc_rrm``) as fixed effects.

        >>> random = RandomState(0)
        >>>
        >>> mean = ones((N, 1))
        >>> pc_rrm = random.randn(N, 4)
        >>> covs = concatenate([mean, pc_rrm], 1)
        >>>
        >>> mtset = MTSet(pheno, F=covs)
        >>> res_null = mtset.fitNull()
        >>>
        >>> print(res_null['conv'][0])
        True
        >>> print('%.2f'%res_null['NLL0'])
        118.36
        >>>
        >>> res = mtset.optimize(snp_set)
        >>>
        >>> print("%.4f" % res['LLR'][0])
        0.1373
        >>> print(res['Cr'])
        [[ 0.0002  0.0019]
         [ 0.0019  0.0207]]
    """

    def __init__(
            self,
            Y=None,
            R=None,
            S_R=None,
            U_R=None,
            traitID=None,
            F=None,
            rank=1):
        from limix_core.gp import GP2KronSum
        from limix_core.gp import GP2KronSumLR
        from limix_core.gp import GP3KronSumLR
        from limix_core.covar import FreeFormCov
        # data
        noneNone = S_R is not None and U_R is not None
        self.bgRE = R is not None or noneNone
        # fixed effect
        msg = 'The current implementation of the full rank mtSet'
        msg += ' does not support covariates.'
        msg += ' We reccommend to regress out covariates and'
        msg += ' subsequently quantile normalize the phenotypes'
        msg += ' to a normal distribution prior to use mtSet.'
        msg += ' This can be done within the LIMIX framework using'
        msg += ' the methods limix.util.preprocess.regressOut and'
        msg += ' limix.util.preprocess.gaussianize'
        assert not (F is not None and self.bgRE), msg
        from limix.util.preprocess import remove_dependent_cols
        if F is not None:
            F = remove_dependent_cols(F)
            A = sp.eye(Y.shape[1])
        else:
            A = None
        # traitID
        if traitID is None:
            traitID = sp.array(['trait %d' % p for p in range(Y.shape[1])])
        self.setTraitID(traitID)
        # init covariance matrices and gp
        Cg = FreeFormCov(Y.shape[1])
        Cn = FreeFormCov(Y.shape[1])
        G = 1. * (sp.rand(Y.shape[0], 1) < 0.2)
        if self.bgRE:
            self._gp = GP3KronSumLR(
                Y=Y,
                Cg=Cg,
                Cn=Cn,
                R=R,
                S_R=S_R,
                U_R=U_R,
                G=G,
                rank=rank)
        else:
            self._gp = GP2KronSumLR(Y=Y, Cn=Cn, G=G, F=F, A=A)
        # null model params
        self.null = None
        # calls itself for column-by-column trait analysis
        self.stSet = None
        self.nullST = None
        self.infoOpt = None
        self.infoOptST = None
        pass

    ##################################################
    # Properties
    ##################################################
    @property
    def N(self):
        return self._gp.covar.dim_r

    @property
    def P(self):
        return self._gp.covar.dim_c

    @property
    def Cr(self):
        return self._gp.covar.Cr

    @property
    def Cg(self):
        return self._gp.covar.Cg

    @property
    def Cn(self):
        return self._gp.covar.Cn

    @property
    def Y(self):
        return self._gp.mean.Y

    @property
    def F(self):
        try:
            return self._gp.mean.F[0]
        except BaseException:
            return None

    @property
    def A(self):
        try:
            return self._gp.mean.A[0]
        except BaseException:
            return None

    @property
    def S_R(self):
        if self.bgRE:
            RV = self._gp.covar.Sr()
        else:
            RV = None
        return RV

    @property
    def U_R(self):
        if self.bgRE:
            RV = self._gp.covar.Lr().T
        else:
            RV = None
        return RV

    #########################################
    # Setters
    #########################################
    @Y.setter
    def Y(self, value):
        assert value.shape[0] == self.N, 'Dimension mismatch'
        assert value.shape[1] == self.P, 'Dimension mismatch'
        self._gp.mean.Y = Y

    def setTraitID(self, traitID):
        self.traitID = traitID

    ################################################
    # Fitting null model
    ###############################################
    def fitNull(
            self,
            cache=False,
            out_dir='./cache',
            fname=None,
            rewrite=False,
            seed=None,
            factr=1e3,
            n_times=10,
            init_method=None,
            verbose=False):
        r"""
        Fit null model

        Args:
            verbose ()
            cache (bool, optional):
                If False (default), the null model is fitted and
                the results are not cached.
                If True, the cache is activated.
                The cache file dir and name can be specified using
                ``hcache`` and ``fname``.
                When ``cache=True``, we distinguish the following cases:

                - if the specified file does not exist,
                  the output of the null model fiting is cached in the file.
                - if the specified file exists and ``rewrite=True``,
                  the cache file is overwritten.
                - if the specified file exists and ``rewrite=False``,
                  the results from the cache file are imported
                  (the null model is not re-fitted).

            out_dir (str, optional):
                output dir of the cache file.
                The default value is "./cache".
            fname (str, optional):
                Name of the cache hdf5 file.
                It must be specified if ``cache=True``.
            rewrite (bool, optional):
                It has effect only if cache `cache=True``.
                In this case, if ``True``,
                the cache file is overwritten in case it exists.
                The default value is ``False``
            factr (float, optional):
                optimization paramenter that determines the accuracy
                of the solution.
                By default it is 1000.
                (see scipy.optimize.fmin_l_bfgs_b for more details).
            verbose (bool, optional):
                verbose flag.

        Returns:
            (dict): dictionary containing:
                - **B** (*ndarray*): estimated effect sizes (null);
                - **Cg** (*ndarray*): estimated relatedness trait covariance (null);
                - **Cn** (*ndarray*): estimated genetic noise covariance (null);
                - **conv** (*bool*): convergence indicator;
                - **NLL0** (*ndarray*): negative loglikelihood (NLL) of the null model;
                - **LMLgrad** (*ndarray*): norm of the gradient of the NLL.
                - **time** (*time*): elapsed time (in seconds).
        """
        from limix_core.gp import GP2KronSum
        from limix_core.gp import GP2KronSumLR
        from limix_core.gp import GP3KronSumLR
        from limix_core.covar import FreeFormCov
        if seed is not None:
            sp.random.seed(seed)

        read_from_file = False
        if cache:
            assert fname is not None, 'MultiTraitSetTest:: specify fname'
            if not os.path.exists(out_dir):
                os.makedirs(out_dir)
            out_file = os.path.join(out_dir, fname)
            read_from_file = os.path.exists(out_file) and not rewrite

        RV = {}
        if read_from_file:
            f = h5py.File(out_file, 'r')
            for key in list(f.keys()):
                RV[key] = f[key][:]
            f.close()
            self.setNull(RV)
        else:
            start = TIME.time()
            if self.bgRE:
                self._gpNull = GP2KronSum(
                    Y=self.Y,
                    F=None,
                    A=None,
                    Cg=self.Cg,
                    Cn=self.Cn,
                    R=None,
                    S_R=self.S_R,
                    U_R=self.U_R)
            else:
                self._gpNull = GP2KronSumLR(
                    self.Y, self.Cn, G=sp.ones(
                        (self.N, 1)), F=self.F, A=self.A)
                # freezes Cg to 0
                n_params = self._gpNull.covar.Cr.getNumberParams()
                self._gpNull.covar.Cr.setParams(1e-9 * sp.ones(n_params))
                self._gpNull.covar.act_Cr = False
            for i in range(n_times):
                params0 = self._initParams(init_method=init_method)
                self._gpNull.setParams(params0)
                conv, info = self._gpNull.optimize(
                    verbose=verbose, factr=factr)
                if conv:
                    break
            if not conv:
                warnings.warn("not converged")
            LMLgrad = (self._gpNull.LML_grad()['covar']**2).mean()
            LML = self._gpNull.LML()
            if self._gpNull.mean.n_terms == 1:
                RV['B'] = self._gpNull.mean.B[0]
            elif self._gpNull.mean.n_terms > 1:
                warning.warn('generalize to more than 1 fixed effect term')
            if self.bgRE:
                RV['params0_g'] = self.Cg.getParams()
            else:
                RV['params0_g'] = sp.zeros_like(self.Cn.getParams())
            RV['params0_n'] = self.Cn.getParams()
            if self.bgRE:
                RV['Cg'] = self.Cg.K()
            else:
                RV['Cg'] = sp.zeros_like(self.Cn.K())
            RV['Cn'] = self.Cn.K()
            RV['conv'] = sp.array([conv])
            RV['time'] = sp.array([TIME.time() - start])
            RV['NLL0'] = sp.array([LML])
            RV['LMLgrad'] = sp.array([LMLgrad])
            RV['nit'] = sp.array([info['nit']])
            RV['funcalls'] = sp.array([info['funcalls']])
            self.null = RV
            from limix.util.util_functions import smartDumpDictHdf5
            if cache:
                f = h5py.File(out_file, 'w')
                smartDumpDictHdf5(RV, f)
                f.close()
        return RV

    def getNull(self):
        """ get null model info """
        return self.null

    def setNull(self, null):
        """ set null model info """
        self.null = null

    ###########################################
    # Fitting alternative model
    ###########################################

    def optimize(
            self,
            G,
            params0=None,
            n_times=10,
            vmax=5,
            factr=1e7,
            verbose=False):
        """
        Fit alternative model

        Args:
            G (ndarray):
                (N, S) genotype values for `N` samples and `S` variants
                (defines the set component)
            params0 (ndarray, optional):
                initial parameter values for optimization.
                By default the parameters are initialized
                by randonmly pertumrning the maximum likelihood
                estimators (MLE) under the null.
            n_times (ndarray, optional):
                maximum number of restarts in case the optimization does not
                converge.
                It has no effect if ``params0`` is provided.
            vmax (float, option):
                maximum value for variance parameters that is accepted.
                If any of the variance MLE is ``>vmax`` then we say that
                the optimization has not converged.
                Default value is 5.
            factr (float, optional):
                optimization paramenter that determines the accuracy
                of the solution.
                By default it is 1e7.
                (see scipy.optimize.fmin_l_bfgs_b for more details).
            verbose (bool, optional):
                verbose flag.

        Returns:
            (dict): dictionary containing:
                - **Cr** (*ndarray*): estimated set trait covariance (alt);
                - **Cg** (*ndarray*): estimated relatedness trait covariance (alt);
                - **Cn** (*ndarray*): estimated genetic noise covariance (alt);
                - **conv** (*ndarray*): convergence indicator;
                - **NLLAlt** (*ndarray*): negative loglikelihood (NLL) of the alternative model;
                - **LMLgrad** (*ndarray*): norm of the gradient of the NLL.
                - **LLR** (*ndarray**): log likelihood ratio statistics;
                - **var** (*ndarray*): (`P`, `K`) ndarray of  variance estimates
                  for the ``P`` traits and the ``K``
                - **time** (*ndarray*): elapsed time (in seconds).
                - **nit** (*ndarray*):
                  number of iterations L-BFGS.
                  (see scipy.optimize.fmin_l_bfgs_b for more details).
        """
        # set params0 from null if params0 is None
        if params0 is None:
            if self.null is None:
                if verbose:
                    print(".. fitting null model")
                self.fitNull()
            if self.bgRE:
                params0 = sp.concatenate(
                    [self.null['params0_g'], self.null['params0_n']])
            else:
                params0 = self.null['params0_n']
            params_was_None = True
        else:
            params_was_None = False
        G *= sp.sqrt(self.N / (G**2).sum())
        self._gp.covar.G = G
        start = TIME.time()
        for i in range(n_times):
            if params_was_None:
                n_params = self.Cr.getNumberParams()
                _params0 = {'covar': sp.concatenate(
                    [1e-3 * sp.randn(n_params), params0])}
            else:
                _params0 = {'covar': params0}
            self._gp.setParams(_params0)
            conv, info = self._gp.optimize(factr=factr, verbose=verbose)
            conv *= self.Cr.K().diagonal().max() < vmax
            conv *= self.getLMLgrad() < 0.1
            if conv or not params_was_None:
                break
        self.infoOpt = info
        if not conv:
            warnings.warn("not converged")
        # return value
        RV = {}
        if self.P > 1:
            RV['Cr'] = self.Cr.K()
            if self.bgRE:
                RV['Cg'] = self.Cg.K()
            RV['Cn'] = self.Cn.K()
        RV['time'] = sp.array([TIME.time() - start])
        RV['params0'] = _params0
        RV['nit'] = sp.array([info['nit']])
        RV['funcalls'] = sp.array([info['funcalls']])
        RV['var'] = self.getVariances()
        RV['conv'] = sp.array([conv])
        RV['NLLAlt'] = sp.array([self.getNLLAlt()])
        RV['LLR'] = sp.array([self.getLLR()])
        RV['LMLgrad'] = sp.array([self.getLMLgrad()])
        return RV

    def getInfoOpt(self):
        """ get information for the optimization """
        return self.infoOpt

    def getVariances(self):
        """
        get variances
        """
        var = []
        var.append(self.Cr.K().diagonal())
        if self.bgRE:
            var.append(self.Cg.K().diagonal())
        var.append(self.Cn.K().diagonal())
        var = sp.array(var)
        return var

    def getNLLAlt(self):
        """
        get negative log likelihood of the alternative
        """
        return self._gp.LML()

    def getLLR(self):
        """
        get log likelihood ratio
        """
        assert self.null is not None, 'null model needs to be fitted!'
        return self.null['NLL0'][0] - self.getNLLAlt()

    def getLMLgrad(self):
        """
        get norm LML gradient
        """
        return (self._gp.LML_grad()['covar']**2).mean()

    def fitNullTraitByTrait(
            self,
            verbose=False,
            cache=False,
            out_dir='./cache',
            fname=None,
            rewrite=False):
        """
        Fit null model trait by trait
        """
        read_from_file = False
        if cache:
            assert fname is not None, 'MultiTraitSetTest:: specify fname'
            if not os.path.exists(out_dir):
                os.makedirs(out_dir)
            out_file = os.path.join(out_dir, fname)
            read_from_file = os.path.exists(out_file) and not rewrite

        RV = {}
        if read_from_file:
            f = h5py.File(out_file, 'r')
            for p in range(self.P):
                trait_id = self.traitID[p]
                g = f[trait_id]
                RV[trait_id] = {}
                for key in list(g.keys()):
                    RV[trait_id][key] = g[key][:]
            f.close()
            self.nullST = RV
        else:
            """ create stSet and fit null column by column returns all info """
            if self.stSet is None:
                y = sp.zeros((self.N, 1))
                self.stSet = MTSet(Y=y, S_R=self.S_R, U_R=self.U_R, F=self.F)
            RV = {}
            for p in range(self.P):
                trait_id = self.traitID[p]
                self.stSet.Y = self.Y[:, p:p + 1]
                RV[trait_id] = self.stSet.fitNull()
            self.nullST = RV
            from limix.util.util_functions import smartDumpDictHdf5
            if cache:
                f = h5py.File(out_file, 'w')
                smartDumpDictHdf5(RV, f)
                f.close()
        return RV

    def optimizeTraitByTrait(self, G, verbose=False, n_times=10, factr=1e3):
        """ Optimize trait by trait """
        assert self.nullST is not None, 'fit null model beforehand'
        RV = {}
        self.infoOptST = {}
        for p in range(self.P):
            trait_id = self.traitID[p]
            self.stSet.Y = self.Y[:, p:p + 1]
            self.stSet.setNull(self.nullST[trait_id])
            RV[trait_id] = self.stSet.optimize(
                G, n_times=n_times, factr=factr, verbose=verbose)
            self.infoOptST[trait_id] = self.stSet.getInfoOpt()
        return RV

    def getInfoOptST(self):
        """ get information on the optimization """
        return self.infoOptST

    def _initParams(self, init_method=None):
        """ internal function for params initialization """
        if self.bgRE:
            if init_method == 'random':
                params0 = {
                    'covar': sp.randn(
                        self._gpNull.covar.getNumberParams())}
            else:
                if self.P == 1:
                    params0 = {'covar': sp.sqrt(0.5) * sp.ones(2)}
                else:
                    cov = 0.5 * sp.cov(self.Y.T) + 1e-4 * sp.eye(self.P)
                    chol = la.cholesky(cov, lower=True)
                    params = chol[sp.tril_indices(self.P)]
                    params0 = {'covar': sp.concatenate([params, params])}
        else:
            if self.P == 1:
                params_cn = sp.array([1.])
            else:
                cov = sp.cov(self.Y.T) + 1e-4 * sp.eye(self.P)
                chol = la.cholesky(cov, lower=True)
                params_cn = chol[sp.tril_indices(self.P)]
            params0 = {'covar': params_cn}
        return params0


if __name__ == '__main__':
    from limix.util.preprocess import covar_rescale

    sp.random.seed(0)

    # generate phenotype
    n = 1000
    p = 4
    f = 10
    Y = sp.randn(n, p)
    X = sp.randn(n, f)
    G = sp.randn(n, f)
    R = sp.dot(X, X.T)
    R = covar_rescale(R)

    mts = mtset(Y, R=R)
    nullMTInfo = mts.fitNull(cache=False)
    mts.optimize(G)
