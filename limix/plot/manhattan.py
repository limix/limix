# Copyright(c) 2014, The LIMIX developers (Christoph Lippert, Paolo Francesco Casale, Oliver Stegle)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import division

from numpy import arange, asarray, cumsum, flipud, log10


def plot_manhattan(df,
                   alpha=None,
                   null_style=dict(alpha=0.1, color='DarkBlue'),
                   alt_style=dict(alpha=0.5, color='Orange'),
                   ax=None):

    import matplotlib.pyplot as plt

    ax = plt.gca() if ax is None else ax

    if 'pos' not in df:
        df['pos'] = arange(df.shape[0])

    df = _abs_pos(df)

    if alpha is None:
        alpha = 0.01 / df.shape[0]

    ytop = -1.2 * log10(min(df['pv'].min(), alpha))

    _plot_chrom_strips(ax, df, ytop)
    _plot_points(ax, df, alpha, null_style, alt_style)
    _set_frame(ax, df, ytop)

    ax.set_ylabel('-log$_{10}$pv')
    ax.set_xlabel('chromosome')

    _set_ticks(ax, _chrom_bounds(df))

    return ax


def _set_frame(ax, df, ytop):
    ax.set_ylim(0, ytop)
    ax.set_xlim(0, df['abs_pos'].max())

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)


def _plot_points(ax, df, alpha, null_style, alt_style):
    ok = df['pv'] >= alpha
    ax.plot(df['abs_pos'][ok], -log10(df['pv'][ok]), '.', ms=5, **null_style)

    ok = df['pv'] < alpha
    ax.plot(df['abs_pos'][ok], -log10(df['pv'][ok]), '.', ms=5, **alt_style)


def _plot_chrom_strips(ax, df, ytop):
    uchroms = df['chrom'].unique()
    for i in range(0, len(uchroms), 2):
        ax.fill_between(
            df['abs_pos'],
            0,
            ytop,
            where=df['chrom'] == uchroms[i],
            facecolor='LightGray',
            linewidth=0,
            alpha=0.5)


def _set_ticks(ax, chrom_bounds):
    n = len(chrom_bounds) - 1
    xticks = asarray([chrom_bounds[i:i + 2].mean() for i in range(n)])
    ax.set_xticks(xticks)
    ax.tick_params(axis='x', which='both', labelsize=6)
    ax.set_xticklabels(arange(1, n + 2))
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')


def _abs_pos(df):
    uchroms = df['chrom'].unique()
    chrom_ends = [df['pos'][df['chrom'] == c].max() for c in uchroms]

    offset = flipud(cumsum(chrom_ends)[:-1])

    df['abs_pos'] = df['pos'].copy()

    uchroms = list(reversed(uchroms))
    for i in range(len(offset)):
        ix = df['chrom'] == uchroms[i]
        df.loc[ix, 'abs_pos'] = df.loc[ix, 'abs_pos'] + offset[i]

    return df


def _chrom_bounds(df):
    uchroms = df['chrom'].unique()
    v = [df['abs_pos'][df['chrom'] == c].min() for c in uchroms]
    return asarray(v + [df['abs_pos'].max()])
