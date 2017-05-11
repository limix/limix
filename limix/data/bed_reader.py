import scipy as sp
import pandas as pd
import copy
from limix.io import read_plink
from sklearn.preprocessing import Imputer


def build_query(idx_start=None,
                idx_end=None,
                pos_start=None,
                pos_end=None,
                chrom=None):
    queries = []

    # gather all queries
    if idx_start is not None:
        query = "i >= %d" % idx_start
        queries.append(query)

    if idx_end is not None:
        query = "i < %d" % idx_end
        queries.append(query)

    if pos_start is not None:
        query = "pos >= %d" % pos_start
        queries.append(query)

    if pos_end is not None:
        query = "pos <% d" % pos_end
        queries.append(query)

    if chrom is not None:
        query = "chrom == '%s'" % str(chrom)
        queries.append(query)

    if len(queries) >= 1:
        query = ' & '.join(queries)
    else:
        query = None

    return query


class BedReader():
    r"""
    Class to read and make queries on plink binary files.

    Parameters
    ----------
    prefix : str
        Path prefix to the set of PLINK files.

    Examples
    --------

    Basic example of usage.

    .. doctest::

        >>> from limix.data.bed_reader import BedReader
        >>> from pandas_plink import example_file_prefix
        >>>
        >>> reader = BedReader(example_file_prefix())
        >>>
        >>> print(reader.getSnpInfo().head())
          chrom         snp   cm    pos a0 a1  i
        0     1  rs10399749  0.0  45162  G  C  0
        1     1   rs2949420  0.0  45257  C  T  1
        2     1   rs2949421  0.0  45413  0  0  2
        3     1   rs2691310  0.0  46844  A  T  3
        4     1   rs4030303  0.0  72434  0  G  4
        >>>
        >>> X, snpinfo = reader.getGenotypes(idx_start=4,
        ...                                  idx_end=10,
        ...                                  pos_start=45200,
        ...                                  pos_end=80000,
        ...                                  chrom=1,
        ...                                  impute=True,
        ...                                  return_snpinfo=True)
        >>>
        >>> print(snpinfo)
          chrom        snp   cm    pos a0 a1  i
        4     1  rs4030303  0.0  72434  0  G  4
        5     1  rs4030300  0.0  72515  0  C  5
        6     1  rs3855952  0.0  77689  G  A  6
        7     1   rs940550  0.0  78032  0  T  7
        >>>
        >>> print(X)
        [[ 2.  2.  2.  2.]
         [ 2.  2.  1.  2.]
         [ 2.  2.  0.  2.]]

    Example of lazt selection using subsampling.

    .. doctest::
        >>> reader_sub = reader.subset_snps(idx_start=4,
        ...                                 idx_end=10,
        ...                                 pos_start=45200,
        ...                                 pos_end=80000,
        ...                                 chrom=1)
        >>> X, snpinfo = reader_sub.getGenotypes(impute=True,
        ...                                      return_snpinfo=True)
        >>> print(X)
        [[ 2.  2.  2.  2.]
         [ 2.  2.  1.  2.]
         [ 2.  2.  0.  2.]]
    """

    def __init__(self, prefix):
        self._prefix = prefix
        self._load()
        self._init_imputer()

    def _load(self):
        (bim, fam, bed) = read_plink(self._prefix, verbose=False)
        self._snpinfo = bim
        self._ind_info = fam
        self._geno = bed

    def _init_imputer(self):
        self._imputer = Imputer(missing_values=3.,
                                strategy='mean',
                                axis=0,
                                copy=False)

    def getSnpInfo(self):
        r"""
        Return pandas dataframe with all variant info.
        """
        return self._snpinfo

    def subset_snps(self,
                    idx_start=None,
                    idx_end=None,
                    pos_start=None,
                    pos_end=None,
                    chrom=None,
                    inplace=False):
        r""" Builds a new bed reader with filtered variants.

        Parameters
        ----------
        idx_start : int, optional
            start idx.
            If not None (default),
            the query 'idx >= idx_start' is considered.
        idx_end : int, optional
            end idx.
            If not None (default),
            the query 'idx < idx_end' is considered.
        pos_start : int, optional
            start chromosomal position.
            If not None (default),
            the query 'pos >= pos_start' is considered.
        pos_end : int, optional
            end chromosomal position.
            If not None (default),
            the query 'pos < pos_end' is considered.
        chrom : int, optional
            chromosome.
            If not None (default),
            the query 'chrom == chrom' is considered.
        copy : bool
            If True, the operation is done in place.
            Default is False.

        Returns
        -------
            R : :class:`limix.BedReader`
                Bed reader with filtered variants.
        """
        # query
        query = build_query(idx_start=idx_start,
                            idx_end=idx_end,
                            pos_start=pos_start,
                            pos_end=pos_end,
                            chrom=chrom)
        X, snpinfo = self._query(query)

        if not inplace:
            # copy (note the first copy is not deep)
            R = copy.copy(self)
            R._ind_info = copy.copy(self._ind_info)
        else:
            # replace
            R = self
        R._geno = X
        R._snpinfo = snpinfo

        return R

    def getGenotypes(self,
                     idx_start=None,
                     idx_end=None,
                     pos_start=None,
                     pos_end=None,
                     chrom=None,
                     impute=False,
                     standardize=False,
                     return_snpinfo=False):
        r""" Query and Load genotype data.

        Parameters
        ----------
        idx_start : int, optional
            start idx.
            If not None (default),
            the query 'idx >= idx_start' is considered.
        idx_end : int, optional
            end idx.
            If not None (default),
            the query 'idx < idx_end' is considered.
        pos_start : int, optional
            start chromosomal position.
            If not None (default),
            the query 'pos >= pos_start' is considered.
        pos_end : int, optional
            end chromosomal position.
            If not None (default),
            the query 'pos < pos_end' is considered.
        chrom : int, optional
            chromosome.
            If not None (default),
            the query 'chrom == chrom' is considered.
        impute : bool, optional
            list of chromosomes.
            If True,
            the missing values in the bed file are mean
            imputed (variant-by-variant).
            If standardize is True, the default value of
            impute is True, otherwise is False.
        standardize : bool, optional
            If True, the genotype values are standardizes.
            The default value is False.
        return_snpinfo : bool, optional
            If True, returns genotype info
            By default is False.

        Returns
        -------
            X : ndarray
                (`N`, `S`) ndarray of queried genotype values
                for `N` individuals and `S` variants.
            snpinfo : :class:`pandas.DataFrame`
                dataframe with genotype info.
                Returned only if ``return_snpinfo=True``.
        """
        if standardize:
            impute = True

        # query
        query = build_query(idx_start=idx_start,
                            idx_end=idx_end,
                            pos_start=pos_start,
                            pos_end=pos_end,
                            chrom=chrom)
        X, snpinfo = self._query(query)

        # compute
        X = X.compute().T

        # impute and standardize
        if impute:
            X = self._imputer.fit_transform(X)

        if standardize:
            X = X.astype(float)
            X -= X.mean(0)
            X /= X.std(0)

        if return_snpinfo:
            return X, snpinfo
        else:
            return X

    def _query(self, query):
        if query is None:
            return self._geno, self._snpinfo
        snpinfo = self._snpinfo.query(query)
        geno = self._geno[snpinfo.i, :]
        return geno, snpinfo
