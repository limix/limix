from pandas import read_csv, MultiIndex

def read_gen(prefix):

    df_sample = read_csv(prefix + '.sample', header=0, sep=' ', skiprows=[1])

    nsamples = df_sample.shape[0]

    col_level0_names = ['snp_id', 'rs_id', 'pos', 'alleleA', 'alleleB']
    col_level1_names = [''] * 5
    for s in df_sample['sample_id']:
        col_level0_names += [s] * 3
        col_level1_names += ['AA', 'AB', 'BB']

    tuples = list(zip(col_level0_names, col_level1_names))
    index = MultiIndex.from_tuples(tuples, names=['first', 'second'])

    df_gen = read_csv(prefix + '.gen', names=index, sep=' ')

    return dict(sample=df_sample, genotype=df_gen)
