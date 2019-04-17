from .._bits.dask import array_shape_reveal


def allele_frequency(X):
    """
    Compute allele frequency from its expectation.

    ``X`` is a matrix whose last dimension correspond to the different set of
    chromosomes. The first dimension represent different variants and the second
    dimension represent the different samples.

    Parameters
    ----------
    X : array_like
        Allele expectations encoded as a variants-by-samples-by-alleles matrix.

    Returns
    -------
    frequenct : ndarray
        Allele frequencies encoded as a variants-by-alleles matrix.
    """
    ploidy = X.shape[-1]
    if X.ndim < 3:
        n = 1
    else:
        n = X.shape[1]
    return X.sum(-2) / (ploidy * n)


def compute_dosage(X, alt=None):
    """
    Compute dosage from allele expectation.

    Parameters
    ----------
    X : array_like
        Allele expectations encoded as a variants-by-samples-by-alleles matrix.
    ref : array_like
        Allele reference of each locus. The allele having the minor allele frequency for
        the provided ``X`` is used as the reference if `None`. Defaults to ``None``.

    Returns
    -------
    dosage : ndarray
        Dosage encoded as a variants-by-samples matrix.

    Example
    -------

    .. doctest::

        >>> from bgen_reader import read_bgen, allele_expectation, example_files
        >>> from bgen_reader import compute_dosage
        >>>
        >>> with example_files("example.32bits.bgen") as filepath:
        ...     bgen = read_bgen(filepath, verbose=False)
        ...     variant_idx = 2
        ...     e = allele_expectation(bgen, variant_idx)
        ...     dosage = compute_dosage(e)
        ...     print(dosage[:5])
        [0.01550294 0.99383543 1.97933958 0.99560547 1.97879027]
        """
    from numpy import asarray

    if alt is None:
        return X[..., -1]
    try:
        return X[alt, :, alt]
    except NotImplementedError:
        alt = asarray(alt, int)
        return asarray(X, float)[alt, :, alt]


def allele_expectation(p, nalleles, ploidy):
    """
    Allele expectation.

    Compute the expectation of each allele from the given probabilities.
    It accepts three shapes of matrices:
    - unidimensional array of probabilities;
    - bidimensional samples-by-alleles probabilities array;
    - and three dimensional variants-by-samples-by-alleles array.

    Parameters
    ----------
    p : array_like
        Allele probabilities.
    nalleles : int
        Number of alleles.
    ploidy : int
        Number of complete sets of chromosomes.

    Returns
    -------
    expectation : ndarray
        Last dimension will contain the expectation of each allele.

    Examples
    --------

    .. doctest::

        >>> from texttable import Texttable
        >>> from bgen_reader import read_bgen, allele_expectation, example_files
        >>>
        >>> sampleid = "sample_005"
        >>> rsid = "RSID_6"
        >>>
        >>> with example_files("example.32bits.bgen") as filepath:
        ...     bgen = read_bgen(filepath, verbose=False)
        ...
        ...     locus = bgen["variants"].query("rsid == '{}'".format(rsid)).index.compute().item()
        ...     sample = bgen["samples"].to_frame().query("id == '{}'".format(sampleid)).index
        ...     sample = sample[0]
        ...
        ...     nalleles = bgen["variants"].loc[locus]["nalleles"]
        ...     ploidy = 2
        ...
        ...     p = bgen["genotype"][locus].compute()["probs"][sample]
        ...     # For unphased genotypes only.
        ...     e = allele_expectation(bgen, locus)[sample]
        ...
        ...     alleles = bgen["variants"].loc[locus]["allele_ids"].compute().item().split(",")
        ...
        ...     tab = Texttable()
        ...
        ...     print(tab.add_rows(
        ...         [
        ...             ["", "AA", "AG", "GG", "E[.]"],
        ...             ["p"] + list(p) + [1.0],
        ...             ["#" + alleles[0], 2, 1, 0, e[0]],
        ...             ["#" + alleles[1], 0, 1, 2, e[1]],
        ...         ]
        ...     ).draw())
        +----+-------+-------+-------+-------+
        |    |  AA   |  AG   |  GG   | E[.]  |
        +====+=======+=======+=======+=======+
        | p  | 0.012 | 0.987 | 0.001 | 1     |
        +----+-------+-------+-------+-------+
        | #A | 2     | 1     | 0     | 1.011 |
        +----+-------+-------+-------+-------+
        | #G | 0     | 1     | 2     | 0.989 |
        +----+-------+-------+-------+-------+
        >>> print("variant: {}".format(rsid))
        variant: RSID_6
        >>> print("sample : {}".format(sampleid))
        sample : sample_005

    Note
    ----
    This function supports unphased genotypes only.
    """
    from numpy import asarray, newaxis

    g = _get_genotypes(ploidy, nalleles)
    c = asarray(_genotypes_to_allele_counts(g), float)
    c = c.T.reshape((1,) * (p.ndim - 1) + (c.shape[1], c.shape[0]))
    p = array_shape_reveal(p)
    return (c * p[..., newaxis, :]).sum(-1)


def _get_genotypes(ploidy, nalleles):
    g = _make_genotypes(ploidy, 1, nalleles)
    g = sorted([list(reversed(i)) for i in g])
    g = [list(reversed(i)) for i in g]
    return g


def _make_genotypes(ploidy, start, end):
    tups = []
    if ploidy == 0:
        return tups
    if ploidy == 1:
        return [[i] for i in range(start, end + 1)]
    for i in range(start, end + 1):
        t = _make_genotypes(ploidy - 1, i, end)
        for ti in t:
            tups += [[i] + ti]
    return tups


def _genotypes_to_allele_counts(genotypes):
    nalleles = genotypes[-1][0]
    counts = []
    for g in genotypes:
        count = [0] * nalleles
        for gi in g:
            count[gi - 1] += 1
        counts.append(count)
    return counts
