Variance decomposition
^^^^^^^^^^^^^^^^^^^^^^

Limix enables flexible fitting of variance component models. Here, we illustrate the
usage of variance component models fit to single genes. It is based on the
`eQTL basics tutorial`_ of limix 1.0, which is now deprecated.

.. _eQTL basics tutorial: https://github.com/limix/limix-tutorials/blob/master/eQTL/eQTL_basics.ipynb

Importing limix
---------------

.. nbplot::

    >>> import limix

The Genetic Model
-----------------

The model considered by the LIMIX variance decomposition module is an extension of the
genetic model employed in standard GWAS, however considers multiple random effects
terms:

.. math::

    \mathbf{y} = \mathbf{F}\boldsymbol{\alpha} + \sum_{x}\mathbf{u}^{(x)} +
            \boldsymbol{\psi},\;\;\;\;
    \mathbf{u}^{(x)}\sim\mathcal{N}
        \left(\mathbf{0},{\sigma^{(x)}}^2\mathbf{R}^{(x)}\right),\;
    \boldsymbol{\Psi}\sim\mathcal{N}\left(\mathbf{0},\sigma_e^2\mathbf{I}_N\right)

where

.. math::

    \begin{eqnarray}
    \mathbf{y}   &=& \text{phenotype vector} \in \mathcal{R}^{N,1} \\
    \mathbf{F}   &=& \text{matrix of $K$ covariates} \in \mathcal{R}^{N,K} \\
    \boldsymbol{\alpha} &=& \text{effect of covariates} \in \mathcal{R}^{K,1} \\
    \mathbf{R}^{(x)}   &=& \text{is the sample covariance matrix for contribution $x$}
                \in \mathcal{R}^{N,N} \\
    \end{eqnarray}

If :math:`\mathbf{R}` is a genetic contribution from set of SNPs :math:`\mathcal{S}`,
with a bit of abuse of notation we can write
:math:`\mathbf{R}= \frac{1}{C}\mathbf{X}_{:,\,\mathcal{S}}{\mathbf{X}_{:,\,\mathcal{S}}}^T`
where
:math:`C=\frac{1}{N}\text{trace}\left(\mathbf{X}_{:,\,\mathcal{S}_i}{\mathbf{X}_{:,\,\mathcal{S}_i}}^T\right)`.

Limix supports an arbitrary number of fixed or random effects to be included in the
model.

Example: cis/trans Variance Decomposition in the glucose condition
------------------------------------------------------------------

Here we use the LIMIX variance decomposition module to quantify the variability in gene
expression explained by proximal (cis) and distal (trans) genetic variation. To do so, we
build a linear mixed model with a fixed effect intercept, two random effects for cis and
trans genetic effects and a noise random effect:

.. math::

    \mathbf{y} = \mathbf{1}_N\mu + \mathbf{u}^{(cis)} + \mathbf{u}^{(trans)}
        + \boldsymbol{\psi},\;\;\;\;
    \mathbf{u}^{(cis)}\sim\mathcal{N}
        \left(\mathbf{0},{\sigma^{(x)}}^2\mathbf{R}^{(cis)}\right), \\
    \mathbf{u}^{(trans)}\sim
        \mathcal{N}\left(\mathbf{0},{\sigma^{(x)}}^2\mathbf{R}^{(trans)}\right),
    \boldsymbol{\Psi}\sim
        \mathcal{N}\left(\mathbf{0},\sigma_e^2\mathbf{I}_N\right)

where :math:`\mathbf{R}^\text{(cis)}` and :math:`\mathbf{R}^\text{(trans)}` are the local
and distal relatedeness matrices, built considering all SNPs in cis and trans (i.e., not
in cis) respectively. As cis region is defined by the 50kb region around each gene.

The gene-model is fitted to gene expression in environment 0 for all genes in the Lysine
Biosynthesis pathway and variance components are averaged thereafter to obtain pathway
based variance components.

.. nbplot::

    >>> from numpy import dot
    >>> import seaborn as sns
    >>> from matplotlib import pyplot as plt
    >>> from matplotlib.ticker import FormatStrFormatter
    >>> from pandas import DataFrame
    >>>
    >>> url = "http://rest.s3for.me/limix/smith08.hdf5.bz2"
    >>> limix.download(url, verbose=False)
    >>> limix.extract("smith08.hdf5.bz2", verbose=False)
    >>> data = limix.io.read_hdf5_limix("smith08.hdf5")
    >>> print(data['phenotype']['col_header'].head())
       environment  gene_ID gene_chrom  gene_end  gene_start gene_strand phenotype_ID  i
    0      0.00000  YOL161C         15     11548       11910           C    YOL161C:0  0
    1      0.00000  YJR107W         10    628319      627333           W    YJR107W:0  1
    2      0.00000  YPL270W         16     32803       30482           W    YPL270W:0  2
    3      0.00000  YGR256W          7   1006108     1004630           W    YGR256W:0  3
    4      0.00000  YDR518W          4   1480153     1478600           W    YDR518W:0  4
    >>> G_all = data["genotype"]["matrix"]
    >>> geno_metadata = data["genotype"]["col_header"]
    >>> K_all = dot(G_all, G_all.T)
    >>> phenotype = data["phenotype"]["matrix"]
    >>> pheno_metadata = data["phenotype"]["col_header"]
    >>>
    >>> # Genes from lysine biosynthesis pathway.
    >>> lysine_group = [
    ...     "YIL094C",
    ...     "YDL182W",
    ...     "YDL131W",
    ...     "YER052C",
    ...     "YBR115C",
    ...     "YDR158W",
    ...     "YNR050C",
    ...     "YJR139C",
    ...     "YIR034C",
    ...     "YGL202W",
    ...     "YDR234W",
    ... ]
    >>> window_size = int(5e5)
    >>>
    >>> res = []

for gene in lysine_group[:2]:
    # Select the row corresponding to gene of interest on environment 0.0.
    query = "(gene_ID == '{}') & (environment == 0.0)".format(gene)
    df = pheno_metadata.query(query)

    pheno_idx = df.i.item()
    gene_pos = df[["gene_chrom", "gene_end", "gene_start"]]
    # Estimated middle point of the gene.
    midpoint = (gene_pos["gene_end"].item() - gene_pos["gene_start"].item()) / 2

    # Window definition.
    start = midpoint - window_size // 2
    end = midpoint + window_size // 2
    geno = geno_metadata.query("(pos >= {}) & (pos <= {})".format(start, end))

    y = phenotype[:, pheno_idx]
    G_cis = G_all[:, geno.i.values]
    K_cis = dot(G_cis, G_cis.T)
    K_trans = limix.qc.normalise_covariance(K_all - K_cis)
    K_cis = limix.qc.normalise_covariance(K_cis)

    # Definition of the model to fit our data from which we extract
    # the relative signal strength.
    glmm = limix.glmm.GLMMComposer(len(y))
    glmm.y = y
    glmm.fixed_effects.append_offset()
    glmm.covariance_matrices.append(K_cis)
    glmm.covariance_matrices.append(K_trans)
    glmm.covariance_matrices.append_iid_noise()
    glmm.fit(verbose=False, progress=False)

    cis_scale = glmm.covariance_matrices[0].scale
    trans_scale = glmm.covariance_matrices[1].scale
    noise_scale = glmm.covariance_matrices[2].scale

    res.append([cis_scale, trans_scale, noise_scale])

res = DataFrame(res, columns=["cis", "trans", "noise"])
res = res.div(res.sum(axis=1), axis=0).mean(axis=0)
res *= 100

ax = sns.barplot(x=res.index, y=res.values)  # doctest: +SKIP
ax.yaxis.set_major_formatter(FormatStrFormatter("%.0f%%"))  # doctest: +SKIP

plt.show()  # doctest: +SKIP


Appendix
--------


This tutorial illustrates the use of limix to anlayse expression datasets. For this
illustration, we consider gene expression levels from a yeast genetics study with freely
available data <cite data-cite="smith2008gene">(Smith & Kruglyak, 2008)</cite>. These
data span 109 individuals with 2,956 marker SNPs and expression levels for 5,493 in
glucose and ethanol growth media respectively.

We start out by discussing how to do QTL mapping, implement models that consider multi
loci and introduce the application of variance component models for single quantitative
traits. Subsequently, these analysis are extended to the corresponding multi-trait
models.
