import pandas as pd
import scipy as sp


def estCumPos(position, offset=0, chrom_len=None, return_chromstart=False):
    r"""
    Compute the cumulative position of variants from position dataframe

    Args:
        position (list or :class:`pandas.DataFrame`):
            positions in chromosome/chromosomal basepair position format.
            It can be specified as

            - list ``[chrom, pos]`` where ``chrom`` and ``pos`` are
              *ndarray* with chromosome values and basepair positions;
            - pandas DataFrame of chromosome values (key='chrom') and
              basepair positions (key='pos').

        chrom_len (ndarray, optional):
            vector with predefined chromosome length.
            By default, the length of the chromosome is taken to be the
            maximum basepair position (key='pos') in ``position`` on that
            chromosome.
        offset (float, optional):
            offset between chromosomes for cumulative position
            (default is 0 bp).
        return_chromstart (bool, optional):
            if ``True``, starting cumulative position of each chromosome
            is also returned (default is ``False``).

    Returns:
        (tuple): tuple containing:
            - **pos_cum** (*ndarray*):
              cumulative positions.
            - **chromstart** (*array_like*):
              starting cumulative positions for each chromosome.
              Returned only if ``return_chromstart=True``.

    Examples
    --------

        This function can be applied on a list of chrom and pos arrays

        .. doctest::

            >>> import scipy as sp
            >>> import pandas as pd
            >>> from limix.util import estCumPos
            >>>
            >>> pos = sp.kron(sp.ones(2), sp.arange(1,5)).astype(int)
            >>> chrom = sp.kron(sp.arange(1,3), sp.ones(4)).astype(int)
            >>>
            >>> pos_cum, chromstart = estCumPos([chrom, pos],
            ...                                 return_chromstart=True)
            >>>
            >>> print(chrom)
            [1 1 1 1 2 2 2 2]
            >>>
            >>> print(pos)
            [1 2 3 4 1 2 3 4]
            >>>
            >>> print(pos_cum)
            [1 2 3 4 5 6 7 8]
            >>>
            >>> print(chromstart)
            [1 5]

        or on a position dataframe:

        .. doctest::

            >>> position = pd.DataFrame(sp.array([chrom, pos]).T,
            ...                         columns=['chrom', 'pos'])
            >>> pos_cum, chromstart = estCumPos(position,
            ...                                 return_chromstart=True)
            >>> position['pos_cum'] = pos_cum
            >>> print(position)
               chrom  pos  pos_cum
            0      1    1        1
            1      1    2        2
            2      1    3        3
            3      1    4        4
            4      2    1        5
            5      2    2        6
            6      2    3        7
            7      2    4        8
    """
    if type(position) != pd.core.frame.DataFrame:

        if type(position) != list:
            raise TypeError('position must be a dataframe or list of arrays')

        chrom = position[0]
        pos = position[1]

        if type(pos) != sp.ndarray or type(chrom) != sp.ndarray:
            raise TypeError('position must be a dataframe or list of arrays')

        if len(pos.shape) != 1 or len(chrom.shape) != 1:
            raise ValueError('pos and chrom should be arrays')

        if len(pos) != len(chrom):
            raise ValueError('pos and chrom should have the same length')

        position = pd.DataFrame(
            sp.array([pos, chrom]).T, columns=['pos', 'chrom'])

    RV = position.copy()

    chromvals = sp.unique(position['chrom'])  # sp.unique is always sorted
    chromstart = sp.zeros_like(chromvals)
    pos_cum = sp.zeros_like(position.shape[0])
    if 'pos_cum' not in position:
        RV["pos_cum"] = sp.zeros_like(position['pos'])
    pos_cum = RV['pos_cum'].values
    to_add = 0
    for i, mychrom in enumerate(chromvals):
        i_chr = position['chrom'] == mychrom
        if chrom_len is None:
            maxpos = position['pos'][i_chr].max() + offset
        else:
            maxpos = chrom_len[i] + offset
        pos_cum[i_chr.values] = to_add + position.loc[i_chr, 'pos']
        chromstart[i] = pos_cum[i_chr.values].min()
        to_add += maxpos

    if return_chromstart:
        return pos_cum, chromstart
    else:
        return pos_cum


def unique_variants(snps, return_idxs=False):
    r"""
    Filters out variants with the same genetic profile.

    Args:
        snps (ndarray):
            (`N`, `S`) ndarray of genotype values for `N` individuals and
            `S` variants.

    Returns:
        ndarray: genotype array with unique variants.

    Examples
    --------

        .. doctest::

            >>> from numpy.random import RandomState
            >>> from numpy import kron, ones
            >>> from limix.util import unique_variants
            >>> from numpy import set_printoptions
            >>> set_printoptions(4)
            >>> random = RandomState(1)
            >>>
            >>> N = 4
            >>> snps = kron(random.randn(N,3)<0., ones((1,2)))
            >>>
            >>> print(snps)
            [[ 0.  0.  1.  1.  1.  1.]
             [ 1.  1.  0.  0.  1.  1.]
             [ 0.  0.  1.  1.  0.  0.]
             [ 1.  1.  0.  0.  1.  1.]]
            >>>
            >>> snps_u = unique_variants(snps)
            >>>
            >>> print(snps_u)
            [[ 0.  1.  1.]
             [ 1.  0.  1.]
             [ 0.  1.  0.]
             [ 1.  0.  1.]]
    """

    _s = sp.dot(sp.rand(snps.shape[0]), snps)

    v, ix = sp.unique(_s, return_index=True)
    idxs_u = sp.sort(ix)

    snps = snps[:, idxs_u]

    if return_idxs:
        return (snps, ix)

    return snps
