from limix.data import build_geno_query


class GIter:
    def __init__(self, geno_reader, start=0, end=None, batch_size=1000):
        self.gr = geno_reader
        self.current = start
        if end is None:
            end = self.gr.getSnpInfo().shape[0]
        self.end = end
        self.batch_size = batch_size

    def __iter__(self):
        return self

    def __next__(self):
        if self.current >= self.end:
            raise StopIteration
        else:
            _end = self.current + self.batch_size
            query = build_geno_query(idx_start=self.current, idx_end=_end)
            rv = self.gr.subset_snps(query)
            self.current = _end
            return rv

    def next(self):  # for python2 compatibility
        return self.__next__()
