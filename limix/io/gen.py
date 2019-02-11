def read(prefix, verbose=True):
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

        >>> from limix.io.gen import read
        >>> from limix import file_example
        >>>
        >>> with file_example(["example.gen", "example.sample"]) as filepaths:
        ...     data = read(filepaths[0][:-4])
        >>>
        >>> print(data['sample'])
               sample_id subject_id  missing  gender  age  age_of_onset  phenotype_1
        sample
        1A0          1A0       W001  0.00000       2    4            -9            0
        1A1          1A1       W002  0.00000       2    4            -9            0
        1A2          1A2       W003  0.00000       2    4            -9            1
        1A3          1A3       W004  0.09000       2    4            -9            1
        1A4          1A4       W005  0.00000       2    4            -9            1
        >>> print(data['genotype'].head())
                  snp_id  rs_id       pos alleleA alleleB 1A0       1A1       1A2       1A3  \
                                                           AA AB BB  AA AB BB  AA AB BB  AA
        candidate
        SA1          SA1  rs001  10000000       A       G   0  0  1   0  0  1   0  0  1   0
        SA2          SA2  rs002  10010000       A       G   0  0  1   0  1  0   1  0  0   0
        SA3          SA3  rs003  10020000       C       T   1  0  0   0  1  0   0  0  1   0
        SA4          SA4  rs004  10030000       G       T   1  0  0   0  1  0   0  0  1   0
        SA5          SA5  rs005  10040000       C       G   0  0  1   0  1  0   1  0  0   0
        <BLANKLINE>
                                    1A4
                        AB       BB  AA       AB       BB
        candidate
        SA1        0.42770  0.57210   0  0.02070  0.97920
        SA2        1.00000  0.00000   1  0.00000  0.00000
        SA3        0.99670  0.00000   0  0.00000  1.00000
        SA4        1.00000  0.00000   0  0.00000  1.00000
        SA5        1.00000  0.00000   1  0.00000  0.00000
    """

    from pandas import read_csv, MultiIndex

    df_sample = read_csv(prefix + ".sample", header=0, sep=" ", skiprows=[1])

    col_level0_names = ["snp_id", "rs_id", "pos", "alleleA", "alleleB"]
    col_level1_names = [""] * 5
    for s in df_sample.iloc[:, 0]:
        col_level0_names += [s] * 3
        col_level1_names += ["AA", "AB", "BB"]

    tuples = list(zip(col_level0_names, col_level1_names))
    index = MultiIndex.from_tuples(tuples, names=["first", "second"])

    df_sample["sample"] = df_sample["sample_id"]
    df_sample = df_sample.set_index("sample")
    df_sample.index.name = "sample"

    df_gen = read_csv(prefix + ".gen", names=index, sep=" ")
    df_gen["candidate"] = df_gen["snp_id"]
    df_gen = df_gen.set_index("candidate")
    df_gen.index.name = "candidate"

    return dict(sample=df_sample, genotype=df_gen)
