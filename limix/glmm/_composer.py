from numpy import all as npall
from numpy import isfinite
from ..likelihood import assert_likelihood_name


class GLMMComposer(object):
    def __init__(self):
        self._liknorm = "normal"
        self._y = None
        self._fixed_effects = FixedEffects()
        self._covariance_matrices = CovarianceMatrices()

    @property
    def likname(self):
        return self._likname

    @likname.setter
    def likname(self, likname):
        assert_likelihood_name(likname)
        self._likname = likname.lower()

    @property
    def y(self):
        return self._y

    @y.setter
    def y(self, y):
        if not npall(isfinite(y)):
            raise ValueError("Phenotype values must be finite.")
        self._y = y

    @property
    def fixed_effects(self):
        return self._fixed_effects

    @property
    def covariance_matrices(self):
        return self._covariance_matrices


class FixedEffects(object):
    def __init__(self):
        self._fixed_effects = []

    def append_offset(self):
        self._fixed_effects.append(None)

    def append(self, m):
        if m.ndim > 2:
            raise ValueError("Fixed-effect has to have between one and two dimensions.")

        if not npall(isfinite(m)):
            raise ValueError("Fixed-effect values must be finite.")

        self._fixed_effects.append(m)


class CovarianceMatrices(object):
    def __init__(self):
        self._covariance_matrices = []

    def append_iid_noise(self):
        self._covariance_matrices.append(None)

    def append(self, K):
        if K.ndim != 2:
            raise ValueError("Covariance-matrix has to have two dimensions.")

        if not npall(isfinite(K)):
            raise ValueError("Covariance-matrix values must be finite.")

        self._covariance_matrices.append(K)

