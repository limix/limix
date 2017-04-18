import scipy as sp
import scipy.linalg as la
import scipy.stats as st

def standardize(Y,in_place=False):
    """
    standardize Y in a way that is robust to missing values
    in_place: create a copy or carry out inplace opreations?
    """
    if in_place:
        YY = Y
    else:
        YY = Y.copy()
    for i in range(YY.shape[1]):
        Iok = ~SP.isnan(YY[:,i])
        Ym = YY[Iok,i].mean()
        YY[:,i]-=Ym
        Ys = YY[Iok,i].std()
        YY[:,i]/=Ys
    return YY

def covar_rescaling_factor(C):
    """
    Returns the rescaling factor for the Gower normalizion on covariance matrix C
    the rescaled covariance matrix has sample variance of 1
    """
    n = C.shape[0]
    P = sp.eye(n) - sp.ones((n,n))/float(n)
    trPCP = sp.trace(sp.dot(P,sp.dot(C,P)))
    r = (n-1) / trPCP
    return r

def covar_rescaling_factor_efficient(C):
    """
    Returns the rescaling factor for the Gower normalizion on covariance matrix C
    the rescaled covariance matrix has sample variance of 1
    """
    n = C.shape[0]
    P = sp.eye(n) - sp.ones((n,n))/float(n)
    CP = C - C.mean(0)[:, sp.newaxis]
    trPCP = sp.sum(P * CP)
    r = (n-1) / trPCP
    return r

def covar_rescale(C):
    """
    Perform Gower normalizion on covariance matrix C
    the rescaled covariance matrix has sample variance of 1
    """
    sf = covar_rescaling_factor(C)
    return sf * C


def toRanks(A):
    """
    converts the columns of A to ranks
    """
    AA=sp.zeros_like(A)
    for i in range(A.shape[1]):
        AA[:,i] = st.rankdata(A[:,i])
    AA=sp.array(sp.around(AA),dtype="int")-1
    return AA

def gaussianize(Y):
    """
    Gaussianize X: [samples x phenotypes]
    - each phentoype is converted to ranks and transformed back to normal using the inverse CDF
    """
    N,P = Y.shape

    YY=toRanks(Y)
    quantiles=(sp.arange(N)+0.5)/N
    gauss = st.norm.isf(quantiles)
    Y_gauss=sp.zeros((N,P))
    for i in range(P):
        Y_gauss[:,i] = gauss[YY[:,i]]
    Y_gauss *= -1
    return Y_gauss

def rankStandardizeNormal(Y):
    """
    Gaussianize X: [samples x phenotypes]
    - each phentoype is converted to ranks and transformed back to normal using the inverse CDF
    """
    return gaussianize(Y)

def regressOut(Y, X, return_b=False):
    """
    regresses out X from Y
    """
    Xd = la.pinv(X)
    b = Xd.dot(Y)
    Y_out = Y-X.dot(b)
    if return_b:
        return Y_out, b
    else:
        return Y_out

def remove_dependent_cols(M, tol=1e-6, display=False):
    """
    Returns a matrix where dependent columsn have been removed
    """
    R = la.qr(M, mode='r')[0][:M.shape[1], :]
    I = (abs(R.diagonal())>tol)
    if sp.any(~I) and display:
        print(('cols ' + str(sp.where(~I)[0]) +
                ' have been removed because linearly dependent on the others'))
        R = M[:,I]
    else:
        R = M.copy()
    return R

def boxcox(X):
    """
    Gaussianize X using the Box-Cox transformation: [samples x phenotypes]

    - each phentoype is brought to a positive schale, by first subtracting the minimum value and adding 1.
    - Then each phenotype transformed by the boxcox transformation
    """
    X_transformed = sp.zeros_like(X)
    maxlog = sp.zeros(X.shape[1])
    for i in range(X.shape[1]):
        i_nan = sp.isnan(X[:,i])
        values = X[~i_nan,i]
        X_transformed[i_nan,i] = X[i_nan,i]
        X_transformed[~i_nan,i], maxlog[i] = st.boxcox(values-values.min()+1.0)
    return X_transformed, maxlog
