Variance decomposition
^^^^^^^^^^^^^^^^^^^^^^

Limix enables flexible fitting of variance component models. Here, we illustrate the
usage of variance component models fit to single genes. It is based on the
`eQTL basics tutorial`_ of limix 1.0, which is now deprecated.

.. _eQTL basics tutorial: https://github.com/limix/limix-tutorials/blob/master/eQTL/eQTL_basics.ipynb

Importing limix
---------------

.. plot::
    :context:

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

.. math::

    \mathbf{R}= \frac{1}{C}
        \mathbf{X}_{:,\,\mathcal{S}}{\mathbf{X}_{:,\,\mathcal{S}}}^T

where

.. math::

    C=\frac{1}{N}
        \text{trace}\left(\mathbf{X}_{:,\,\mathcal{S}_i}
        {\mathbf{X}_{:,\,\mathcal{S}_i}}^T\right).

Limix supports an arbitrary number of fixed or random effects to be included in the
model.

Example: cis/trans Variance Decomposition in the glucose condition
------------------------------------------------------------------

Here we use the LIMIX variance decomposition module to quantify the variability in gene
expression explained by proximal (cis) and distal (trans) genetic variation. To do so,
we build a linear mixed model with a fixed effect intercept, two random effects for cis
and trans genetic effects and a noise random effect:

.. math::

    \mathbf{y} = \mathbf{1}_N\mu + \mathbf{u}^{(cis)} + \mathbf{u}^{(trans)}
        + \boldsymbol{\psi},\;\;\;\;
    \mathbf{u}^{(cis)}\sim\mathcal{N}
        \left(\mathbf{0},{\sigma^{(x)}}^2\mathbf{R}^{(cis)}\right), \\
    \mathbf{u}^{(trans)}\sim
        \mathcal{N}\left(\mathbf{0},{\sigma^{(x)}}^2\mathbf{R}^{(trans)}\right),
    \boldsymbol{\Psi}\sim
        \mathcal{N}\left(\mathbf{0},\sigma_e^2\mathbf{I}_N\right)

where :math:`\mathbf{R}^\text{(cis)}` and :math:`\mathbf{R}^\text{(trans)}` are the
local and distal relatedeness matrices, built considering all SNPs in cis and trans
(i.e., not in cis) respectively. As cis region is defined by the 50kb region around
each gene.

The gene-model is fitted to gene expression in environment 0 for all genes in the Lysine
Biosynthesis pathway and variance components are averaged thereafter to obtain pathway
based variance components.

.. plot::
    :context:

    >>> from numpy import dot
    >>> import seaborn as sns
    >>> from matplotlib.ticker import FormatStrFormatter
    >>> from pandas import DataFrame
    >>>
    >>> url = "http://rest.s3for.me/limix/smith08.hdf5.bz2"
    >>> limix.sh.download(url, verbose=False)
    >>> filename = limix.sh.extract("smith08.hdf5.bz2", verbose=False)
    >>> data = limix.io.hdf5.read_limix(filename)
    >>> Y = data['phenotype']
    >>> G_all = data['genotype']
    >>> print(Y)
    <xarray.DataArray 'phenotype' (sample: 109, outcome: 10986)>
    array([[-0.037339, -0.078165,  0.042936, ...,  0.095596, -0.132385, -0.274954],
           [-0.301376,  0.066055,  0.338624, ..., -0.142661, -0.238349,  0.732752],
           [ 0.002661,  0.121835, -0.137064, ..., -0.144404,  0.257615,  0.015046],
           ...,
           [-0.287339,  0.351835,  0.072936, ...,  0.097339, -0.038349,  0.162752],
           [-0.577339,  0.011835, -0.007064, ...,  0.135596,  0.107615,  0.245046],
           [-0.277339,  0.061835,  0.132936, ...,  0.015596, -0.142385, -0.124954]])
    Coordinates:
      * sample        (sample) int64 0 1 2 3 4 5 6 7 ... 102 103 104 105 106 107 108
        environment   (outcome) float64 0.0 0.0 0.0 0.0 0.0 ... 1.0 1.0 1.0 1.0 1.0
        gene_ID       (outcome) object 'YOL161C' 'YJR107W' ... 'YLR118C' 'YBR242W'
        gene_chrom    (outcome) object '15' '10' '16' '7' '4' ... '3' '10' '12' '2'
        gene_end      (outcome) int64 11548 628319 32803 ... 315049 384726 705381
        gene_start    (outcome) int64 11910 627333 30482 ... 315552 385409 704665
        gene_strand   (outcome) object 'C' 'W' 'W' 'W' 'W' ... 'W' 'W' 'C' 'C' 'W'
        phenotype_ID  (outcome) object 'YOL161C:0' 'YJR107W:0' ... 'YBR242W:1'
    Dimensions without coordinates: outcome
    >>> K_all = dot(G_all, G_all.T)
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
    >>>
    >>> for gene in lysine_group[:2]:
    ...     # Select the row corresponding to gene of interest on environment 0.0.
    ...     query = "(gene_ID == '{}') & (environment == 0.0)".format(gene)
    ...     df = Y[:, (Y["gene_ID"] == gene) & (Y["environment"] == 0.0)]
    ...
    ...     # Estimated middle point of the gene.
    ...     midpoint = (df["gene_end"].item() - df["gene_start"].item()) / 2
    ...
    ...     # Window definition.
    ...     start = midpoint - window_size // 2
    ...     end = midpoint + window_size // 2
    ...     geno = G_all[:, (G_all["pos"] >= start) & (G_all["pos"] <= end)]
    ...
    ...     y = df
    ...     G_cis = G_all[:, geno.candidate]
    ...     K_cis = dot(G_cis, G_cis.T)
    ...     K_trans = limix.qc.normalise_covariance(K_all - K_cis)
    ...     K_cis = limix.qc.normalise_covariance(K_cis)
    ...
    ...     # Definition of the model to fit our data from which we extract
    ...     # the relative signal strength.
    ...     glmm = limix.glmm.GLMMComposer(len(y))
    ...     glmm.y = y
    ...     glmm.fixed_effects.append_offset()
    ...     glmm.covariance_matrices.append(K_cis)
    ...     glmm.covariance_matrices.append(K_trans)
    ...     glmm.covariance_matrices.append_iid_noise()
    ...     glmm.fit(verbose=False)
    ...
    ...     cis_scale = glmm.covariance_matrices[0].scale
    ...     trans_scale = glmm.covariance_matrices[1].scale
    ...     noise_scale = glmm.covariance_matrices[2].scale
    ...
    ...     res.append([cis_scale, trans_scale, noise_scale])
    >>>
    >>>
    >>> res = DataFrame(res, columns=["cis", "trans", "noise"])
    >>> res = res.div(res.sum(axis=1), axis=0).mean(axis=0)
    >>> res *= 100
    >>>
    >>> ax = sns.barplot(x=res.index, y=res.values)
    >>> ax.yaxis.set_major_formatter(FormatStrFormatter("%.0f%%"))
    >>>
    >>> limix.plot.show()


We then remove the temporary files.

.. plot::
    :context:

    >>> limix.sh.remove("smith08.hdf5.bz2")
    >>> limix.sh.remove("smith08.hdf5")
