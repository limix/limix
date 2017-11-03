import limix
from pandas import DataFrame
from numpy.random import RandomState

random = RandomState(1)

pv0 = random.rand(10000)
pv1 = random.rand(10000)

pv = list(pv0) + list(pv1)
label = ['label0'] * len(pv0) + ['label1'] * len(pv1)
df = DataFrame(data=dict(pv=pv, label=label))

limix.plot.qqplot(df)
limix.plot.show()
