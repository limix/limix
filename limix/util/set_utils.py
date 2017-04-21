from pandas_plink import read_plink
from optparse import OptionParser
import pandas as pd
import numpy as np


def sets_from_bed(
        bim,
        size=50000,
        step=None,
        chrom=None,
        minSnps=1,
        maxSnps=None):
    r"""
    Builds variant-sets from a bim dataframe considering
    a sliding window approach.

    Args:
        bim (:class:`pandas.DataFrame`):
            bim dataframe from ``:func:pandas_plink.read_plink``.
        size (int, optional):
            set size (in base pairs).
            The default value is 50000.
        step (int, optional):
            sliding-window step.
            The default value is ``size/2``.
        chrom (int, optional):
            can be set to extract set file for
            an only chromosome.
            By default a genome-wide set file
            is built.
        minSnps (int, optional):
            minimin number of variants.
            Sets with number of variants that is
            lower than ``minSnps`` are excluded.
            Default value is 1.
        maxSnps (int, optional):
            maximum number of variants.
            Default value is ``numpy.inf``.

    Returns:
        :class:`pandas.DataFrame`:
            dataframe of variant-sets. It has columns:

            - `"chr"`: chromosome
            - `"start"`: start position
            - `"end"`: end position
            - `"nsnps"`: number of variants in the region
    """
    if step is None:
        step = int(0.5 * size)

    if maxSnps is None:
        maxSnps = np.inf

    if chrom is None:
        chroms = pd.unique(bim['chrom'])
    else:
        chroms = [chrom]

    out = []
    for _c in chroms:

        # calc start and end of regions
        pos = bim.query("chrom=='%s'" % _c)['pos'].values
        start = np.arange(pos.min(), pos.max(), step)
        end = start + size

        # filter based on minSnps
        nsnps = calc_nsnps(pos, start, end)
        Iok = (nsnps >= minSnps) & (nsnps <= maxSnps)
        start = start[Iok]
        end = end[Iok]
        nsnps = nsnps[Iok]

        # build array and append
        chrom = np.repeat(_c, start.shape[0])
        _out = np.array([chrom, start, end, nsnps], dtype=object).T
        out.append(_out)

    # concatenate and export
    out = np.concatenate(out)

    # convert to pandas dataframe and export
    out = pd.DataFrame(out, columns=['chr', 'start', 'end', 'nsnps'])
    return out


def annotate_sets(sets, bim, minSnps=1, maxSnps=None):
    r"""
    Helper function to annotate and filter variant-sets.

    Provided a list of variant-sets and the bim of the bed file to
    analyse, it computes the number of variants within each set
    and filters them.

    Args:
        sets (:class:`pandas.DataFrame`):
            dataframe defining the variant-sets.
            It should contain the columns:

            - `"chr"`: chromosome
            - `"start"`: start position
            - `"end"`: end position

        bim (:class:`pandas.DataFrame`):
            bim dataframe from ``:func:pandas_plink.read_plink``.
        minSnps (int, optional):
            minimin number of variants.
            Sets with number of variants that is
            lower than ``minSnps`` are excluded.
            Default value is 1.
        maxSnps (int, optional):
            maximum number of variants.
            Default value is ``numpy.inf``.

    Returns:
        :class:`pandas.DataFrame`:
            dataframe of variant-sets. It has columns:

            - `"chr"`: chromosome
            - `"start"`: start position
            - `"end"`: end position
            - `"nsnps"`: number of variants in the region
    """
    if maxSnps is None:
        maxSnps = np.inf

    out = []
    chrom = sets['chr'].values
    uchroms = pd.unique(chrom)
    for _c in uchroms:

        # filter on chromosome
        pos = bim.query("chrom=='%s'" % _c)['pos'].values
        Ichr = chrom == _c
        start = sets['start'].values[Ichr]
        end = sets['end'].values[Ichr]

        # filter based on minSnps
        nsnps = calc_nsnps(pos, start, end)
        Iok = (nsnps >= minSnps) & (nsnps <= maxSnps)
        start = start[Iok]
        end = end[Iok]
        nsnps = nsnps[Iok]

        # build array and append
        chrom = np.repeat(_c, start.shape[0])
        _out = np.array([chrom, start, end, nsnps], dtype=object).T
        out.append(_out)

    # concatenate and export
    out = np.concatenate(out)

    # convert to pandas dataframe and export
    out = pd.DataFrame(out, columns=['chr', 'start', 'end', 'nsnps'])
    return out


def calc_nsnps(pos, start, end):
    nsnps = np.zeros(start.shape[0], dtype=int)
    for i, r in enumerate(zip(start, end)):
        nsnps[i] = ((pos >= r[0]) & (pos < r[1])).sum()
    return nsnps


if __name__ == '__main__':

    import os

    bedpath = 'data/chrom22_subsample20_maf0.10'

    if not os.path.exists(bedpath + '.bim'):
        os.system('wget http://www.ebi.ac.uk/~casale/data.zip')
        os.system('unzip data.zip')

    (bim, fam, G) = read_plink(bedpath)

    # sets_from_bed
    sets = sets_from_bed(bim)

    # annotate_sets
    sets1 = annotate_sets(sets, bim, minSnps=100, maxSnps=200)
