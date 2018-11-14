class BaseMean(object):
    def __init__(self, mean):
        self._mean = mean

    @property
    def name(self):
        return self._mean.name

    @name.setter
    def name(self, v):
        self._mean.name = v

    def value(self):
        return self._mean.feed().value()

    def gradient(self):
        return self._mean.feed().gradient()

    def variables(self):
        return self._variables()

    def __str__(self):
        return str(self._mean)


class OffsetMean(BaseMean):
    def __init__(self, mean):
        BaseMean.__init__(self, mean)

    @property
    def offset(self):
        return self._mean.offset

    @offset.setter
    def offset(self, v):
        self._mean.offset = v


class LinearMean(BaseMean):
    def __init__(self, mean):
        BaseMean.__init__(self, mean)

    @property
    def effsizes(self):
        return self._mean.effsizes

    @effsizes.setter
    def effsizes(self, v):
        from numpy import ascontiguousarray

        self._mean.effsizes = ascontiguousarray(v)


class SumMean(BaseMean):
    def __init__(self, mean):
        BaseMean.__init__(self, mean)

    def value(self):
        return self._mean.feed().value()
