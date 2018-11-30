class BaseCov(object):
    def __init__(self, cov):
        self._cov = cov

    @property
    def name(self):
        return self._cov.name

    @name.setter
    def name(self, v):
        self._cov.name = v

    def value(self):
        return self._cov.feed().value()

    def gradient(self):
        return self._cov.feed().gradient()

    def variables(self):
        return self._variables()

    def __str__(self):
        return str(self._cov)


class EyeCov(BaseCov):
    def __init__(self, cov):
        BaseCov.__init__(self, cov)

    @property
    def scale(self):
        return float(self._cov.scale)

    @scale.setter
    def scale(self, v):
        self._cov.scale = v


class GivenCov(BaseCov):
    def __init__(self, cov):
        BaseCov.__init__(self, cov)

    @property
    def scale(self):
        return float(self._cov.scale)

    @scale.setter
    def scale(self, v):
        self._cov.scale = v


class SumCov(BaseCov):
    def __init__(self, cov):
        BaseCov.__init__(self, cov)

    def value(self):
        return self._cov.feed().value()
