**********************
Variance decomposition
**********************

The genetic model
^^^^^^^^^^^^^^^^^

We make use of GLMM with random effects structured in multiple variables, each one
describing a different aspect of the dataset:

.. math::

    \mathbf y = \mathbf M\boldsymbol\beta
        + \sum_j\mathbf X^{(j)}\mathbf u^{(j)} + \boldsymbol\epsilon,

.. math::

    \mathbf u^{(j)} \sim \mathcal N(\mathbf 0, v_j\mathbf I_{j}), ~~\text{and}~~
    \boldsymbol\epsilon\sim\mathcal N(\mathbf 0, v_{j+1}\mathbf I_{j+1}),

where

.. math::

    \begin{eqnarray}
    \mathbf y   &=& \text{phenotype} \in \mathcal R^{n\times 1} \\
    \mathbf M   &=& \text{matrix of covariates} \in \mathcal R^{n\times m} \\
    \boldsymbol\beta &=& \text{covariate effect-sizes} \in \mathcal R^{m\times 1} \\
    \mathbf X^{(j)}   &=& \text{matrix of attributes} \in \mathcal R^{n\times p_j}
    \end{eqnarray}

Example: glucose condition
^^^^^^^^^^^^^^^^^^^^^^^^^^

Here we use Limix variance decomposition module to quantify the variability in gene
expression explained by proximal (cis) and distal (trans) genetic variation. To do so,
we build a linear mixed model with a intercept, two random effects for cis
and trans genetic effects, and a noise random effect.

Lets first download the dataset.

.. plot::
    :context: reset
    :include-source:

    >>> import limix
    >>>
    >>> url = "http://rest.s3for.me/limix/smith08.hdf5.bz2"
    >>> limix.sh.download(url, verbose=False)
    >>> filename = limix.sh.extract("smith08.hdf5.bz2", verbose=False)
    >>> # This dataset in the old limix format.
    >>> data = limix.io.hdf5.read_limix(filename)
    >>> Y = data['phenotype']
    >>> G_all = data['genotype']

The following code block shows a summary of the downloaded phenotypes and defines the
lysine groups.

.. note::

    The phenotype variable ``Y`` is of type :class:`xarray.DataArray`. ``Y`` has
    two dimensions and multiple coordinates associated with them.

.. plot::
    :context:
    :include-source:

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

We will compute the relationship matrix ``K_all`` considering all SNPs and define
the cis region size ``window_size`` in base pairs.
Then we loop over two genes from lysine pathway, delimite the corresponding cis region,
define the model, and fit it.

.. plot::
    :context:
    :include-source:

    >>> from numpy import dot
    >>>
    >>> K_all = dot(G_all, G_all.T)
    >>> window_size = int(5e5)
    >>>
    >>> variances = []
    >>>
    >>> # We loop over the first two groups only.
    >>> for gene in lysine_group[:2]:
    ...
    ...     # Select the row corresponding to gene of interest on environment 0.0.
    ...     y = Y[:, (Y["gene_ID"] == gene) & (Y["environment"] == 0.0)]
    ...
    ...     # Estimated middle point of the gene.
    ...     midpoint = (y["gene_end"].item() - y["gene_start"].item()) / 2
    ...
    ...     # Window definition.
    ...     start = midpoint - window_size // 2
    ...     end = midpoint + window_size // 2
    ...     geno = G_all[:, (G_all["pos"] >= start) & (G_all["pos"] <= end)]
    ...
    ...     G_cis = G_all[:, geno.candidate]
    ...     K_cis = dot(G_cis, G_cis.T)
    ...     # Normalising the covariances is important for comparing their relative
    ...     # overall variances.
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
    ...     variances.append([cis_scale, trans_scale, noise_scale])

We now plot the results.

.. plot::
    :context:
    :include-source:

    >>> import seaborn as sns
    >>> from matplotlib.ticker import FormatStrFormatter
    >>> from pandas import DataFrame
    >>>
    >>> variances = DataFrame(variances, columns=["cis", "trans", "noise"])
    >>> variances = variances.div(variances.sum(axis=1), axis=0).mean(axis=0)
    >>> variances = variances * 100
    >>>
    >>> ax = sns.barplot(x=variances.index, y=variances.values)
    >>> ax.yaxis.set_major_formatter(FormatStrFormatter("%.0f%%"))
    >>>
    >>> limix.plot.show()

And remove temporary files.

.. plot::
    :context: close-figs
    :include-source:

    >>> limix.sh.remove("smith08.hdf5.bz2")
    >>> limix.sh.remove("smith08.hdf5")
