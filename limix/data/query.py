import scipy as sp

class QueryList():
    r"""
    This class helps build str queries for pandas dataframes.

    Parameters
    ----------
    queries : list or str
        list of string queries or single string query

    Examples
    --------

    .. doctest::

        >>> from limix.data import QueryList
        >>>
        >>> ql = QueryList()
        >>> ql.append('i >= 4')
        >>> ql.append('i < 10')
        >>> ql.append('pos >= 45200')
        >>> ql.append('pos < 80000')
        >>> ql.append('chrom == 1')
        >>> print(ql.make_str())
        i >= 4 & i < 10 & pos >= 45200 & pos < 80000 & chrom == 1
    """

    def __init__(self, queries=None):
        self.queries = []
        if queries is not None:
            if type(queries) in [list, str]:
                raise TypeError('queries must be a list or a str')

            if type(queries)==list:
                for query in queries:
                    self.append(query)
            else:
                self.append(query)

    def append(self, query):
        r"""
        Append a new query

        Parameters
        ----------
        query : str
            query to append
        """
        if type(query)!=str:
            raise TypeError('single queries must be strings')
        self.queries.append(query)

    def make_str(self):
        r"""
        Produce a str query joining all queries

        Returns
        -------
            str_query : str 
        """
        if len(self.queries) >= 1:
            str_query = ' & '.join(self.queries)
        else:
            str_query = None
        return str_query
