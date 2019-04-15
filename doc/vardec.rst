**********************
Variance decomposition
**********************

Single-trait decomposition 
==========================

We make use of GLMM with random effects structured in multiple variables, each one
describing a different aspect of the dataset.  For a normally distributed phenotype, we
use the model

.. math::

    𝐲 = 𝐌𝛃 + ∑ⱼ𝐮ⱼ, ~~\text{where}~~ 𝐮ⱼ ∼ 𝓝(𝟎, 𝓋ⱼ𝙺ⱼ).

For non-normally distributed phenotype, the model is given by


.. math::

    𝐳 = 𝐌𝛃 + ∑ⱼ𝐮ⱼ, ~~\text{where}~~~~~~~~~~~~~~~~~~\\
    𝐮ⱼ ∼ 𝓝(𝟎, 𝓋ⱼ𝙺ⱼ) ~~\text{and}~~ yᵢ|𝐳 ∼ 𝙴𝚡𝚙𝙵𝚊𝚖(𝜇ᵢ=g(zᵢ)).

The parameters 𝛃 and 𝓋ⱼ are fit via maximum likelihood.

Example: glucose condition
--------------------------

Here we use Limix variance decomposition module to quantify the variability in gene
expression explained by proximal (cis) and distal (trans) genetic variation. To do so,
we build a linear mixed model with an intercept, two random effects for cis and trans
genetic effects, and a noise random effect.

Lets first download the dataset.

.. plot::
    :context: reset
    :include-source:

    >>> import limix
    >>>
    >>> url = "http://rest.s3for.me/limix/smith08.hdf5.bz2"
    >>> filepath = limix.sh.download(url, verbose=False)
    >>> filepath = limix.sh.extract(filepath, verbose=False)
    >>> # This dataset in the old limix format.
    >>> data = limix.io.hdf5.read_limix(filepath)
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

We will compute the relationship matrix ``K_all`` considering all SNPs and define the
cis region size ``window_size`` in base pairs.  Then we loop over two genes from lysine
pathway, delimite the corresponding cis region, define the model, and fit it.

.. plot::
    :context:
    :include-source:

    >>> from numpy import dot
    >>>
    >>> K_all = dot(G_all, G_all.T)
    >>> window_size = int(5e5)
    >>>
    >>> vardecs = []
    >>>
    >>> # We loop over the first two groups only.
    >>> for gene in lysine_group[:2]:
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
    ...     K_trans = K_all - K_cis
    ...
    ...     # Definition of the model to fit our data from which we extract
    ...     # the relative signal strength.
    ...     vardec = limix.vardec.VarDec(y, "normal")
    ...     vardec.append(K_cis, "cis")
    ...     vardec.append(K_trans, "trans")
    ...     vardec.append_iid("noise")
    ...     vardec.fit(verbose=False)
    ...     vardecs.append(vardec)

We show a summary of each decomposition.

.. plot::
    :context:
    :include-source:

    >>> print(vardecs[0])
    Variance decomposition
    ======================
    <BLANKLINE>
    𝐲 ~ 𝓝(𝙼𝜶, 0.018⋅𝙺 + 0.047⋅𝙺 + 0.066⋅𝙸)
    >>> print(vardecs[1])
    Variance decomposition
    ======================
    <BLANKLINE>
    𝐲 ~ 𝓝(𝙼𝜶, 0.197⋅𝙺 + 0.087⋅𝙺 + 0.149⋅𝙸)

We now plot the results.

.. plot::
    :context:
    :include-source:

    >>> vardecs[0].plot()

.. plot::
    :context:
    :include-source:

    >>> vardecs[1].plot()

And remove temporary files.

.. plot::
    :context: close-figs
    :include-source:

    >>> limix.sh.remove("smith08.hdf5.bz2")
    >>> limix.sh.remove("smith08.hdf5")
