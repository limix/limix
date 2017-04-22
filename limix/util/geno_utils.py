def estCumPos(position, offset=0, chrom_len=None, return_chromstart=False):
    r"""
    Helper function to compute the cumulative position of variants.

    Args:
        position (:class:`pandas.DataFrame`):
            pandas DataFrame of basepair positions (key='pos') and
            chromosome values (key='chrom').
            The DataFrame will be updated with field 'pos_cum'.
        chrom_len (ndarray, optional):
            vector with predefined chromosome length
            By default, the length of the chromosome is assumed to be the
            maximum basepair position (key='pos') in ``position`` on that
            chromosome.
        offset (float, optional):
            offset between chromosomes for cumulative position (default 0 bp).
        return_chromstart (bool, optional):
            if True, starting cumulative position of each chromosome is also
            returned.

    Returns:
        (tuple): tuple containing:
            - **position** (*:class:`pandas.DataFrame`*):
              augmented position dataframe with cumulative positions;
              contains null test statistcs of mtSet, iSet, and iSet-GxC tests.
            - **chrom_poscum** (*array_like*):
              starting cumulative positions for each chromosome.

    Examples
    --------
        .. doctest::

            >>> import scipy as sp
            >>> import pandas as pd
            >>> from limix.util import estCumPos
            >>> 
            >>> pos = sp.kron(sp.ones(2), sp.arange(1,5)).astype(int)
            >>> chrom = sp.kron(sp.arange(1,3), sp.ones(4)).astype(int)
            >>> position = pd.DataFrame(sp.array([pos, chrom]).T, columns=['pos', 'chrom'])
            >>> print(position)
               pos  chrom
            0    1      1
            1    2      1
            2    3      1
            3    4      1
            4    1      2
            5    2      2
            6    3      2
            7    4      2
            >>>
            >>> position, chromStart = estCumPos(position, return_chromstart=True)
            >>> print(position)
               pos  chrom  pos_cum
            0    1      1        1
            1    2      1        2
            2    3      1        3
            3    4      1        4
            4    1      2        5
            5    2      2        6
            6    3      2        7
            7    4      2        8
            >>>
            >>> print(chromStart)
            [1 5]
    """

    RV = position.copy()
    chromvals =  sp.unique(position['chrom'])# sp.unique is always sorted
    chrom_poscum= sp.zeros_like(chromvals)#get the starting position of each Chrom
    pos_cum= sp.zeros_like(position.shape[0])
    if not 'pos_cum' in position:
        RV["pos_cum"]= sp.zeros_like(position['pos'])#get the cum_pos of each variant.
    pos_cum=RV['pos_cum'].values
    to_add=0
    for i,mychrom in enumerate(chromvals):
        i_chr=position['chrom']==mychrom
        if chrom_len is None:
            maxpos = position['pos'][i_chr].max()+offset
        else:
            maxpos = chrom_len[i]+offset
        pos_cum[i_chr.values]=to_add+position.loc[i_chr,'pos']
        chrom_poscum[i] = pos_cum[i_chr.values].min()
        to_add+=maxpos      

    if return_chromstart:
        return RV, chrom_poscum
    else:
        return RV

