def read_gen(prefix):
    r"""
    Read GEN files into Pandas data frames.

    Parameters
    ----------
    prefix : str
        Path prefix to the set of GEN files.

    Returns
    -------
    sample : dask dataframe
    genotype : dask dataframe

    Examples
    --------
    .. doctest::

        >>> from limix.io import read_gen
        >>> from limix.io.examples import gen_file_example
        >>>
        >>> data = read_gen(gen_file_example())
        >>> print(data['sample']) #doctest: +NORMALIZE_WHITESPACE
          sample_id subject_id  missing  gender  age  age_of_onset  phenotype_1
        0       1A0       W001     0.00       2    4            -9            0
        1       1A1       W002     0.00       2    4            -9            0
        2       1A2       W003     0.00       2    4            -9            1
        3       1A3       W004     0.09       2    4            -9            1
        4       1A4       W005     0.00       2    4            -9            1
        >>> print(data['genotype'].head()) #doctest: +NORMALIZE_WHITESPACE
                snp_id  rs_id       pos alleleA alleleB 1A0       1A1       1A2       1A3                 1A4
                                                AA AB BB  AA AB BB  AA AB BB  AA      AB      BB  AA      AB      BB
        0    SA1  rs001  10000000       A       G   0  0  1   0  0  1   0  0  1   0  0.4277  0.5721   0  0.0207  0.9792
        1    SA2  rs002  10010000       A       G   0  0  1   0  1  0   1  0  0   0  1.0000  0.0000   1  0.0000  0.0000
        2    SA3  rs003  10020000       C       T   1  0  0   0  1  0   0  0  1   0  0.9967  0.0000   0  0.0000  1.0000
        3    SA4  rs004  10030000       G       T   1  0  0   0  1  0   0  0  1   0  1.0000  0.0000   0  0.0000  1.0000
        4    SA5  rs005  10040000       C       G   0  0  1   0  1  0   1  0  0   0  1.0000  0.0000   1  0.0000  0.0000
    """

    from pandas import read_csv, MultiIndex

    df_sample = read_csv(prefix + ".sample", header=0, sep=" ", skiprows=[1])

    nsamples = df_sample.shape[0]

    col_level0_names = ["snp_id", "rs_id", "pos", "alleleA", "alleleB"]
    col_level1_names = [""] * 5
    for s in df_sample["sample_id"]:
        col_level0_names += [s] * 3
        col_level1_names += ["AA", "AB", "BB"]

    tuples = list(zip(col_level0_names, col_level1_names))
    index = MultiIndex.from_tuples(tuples, names=["first", "second"])

    df_gen = read_csv(prefix + ".gen", names=index, sep=" ")

    return dict(sample=df_sample, genotype=df_gen)
