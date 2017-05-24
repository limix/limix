def query_and(*queries):
    """
    Given multiple queries it joins them using the & operator.

    Examples
    --------
    .. doctest::
        >>> from limix.data import query_and
        >>>
        >>> query1 = 'a < b'
        >>> query2 = 'b >= c'
        >>> query3 = 'c == 1'
        >>> query = query_and(query1, query2, query3)
        >>>
        >>> print(query)
        a < b & b >= c & c == 1
    """
    if len(queries) >= 1:
        query = ' & '.join(queries)
    else:
        query = None
    return query


def build_geno_query(idx_start=None,
                     idx_end=None,
                     pos_start=None,
                     pos_end=None,
                     chrom=None):
    """
    helper function to build genotype queries.

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

    Returns
    -------
        query : str

    Examples
    --------
    .. doctest::
        >>> from limix.data import build_geno_query
        >>>
        >>> query = build_geno_query(idx_start=4,
        ...                          idx_end=10,
        ...                          pos_start=45200,
        ...                          pos_end=80000,
        ...                          chrom=1)
        >>>
        >>> print(query)
        i >= 4 & i < 10 & pos >= 45200 & pos < 80000 & chrom == '1'
    """
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

    return query_and(*queries)
