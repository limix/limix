******************
Statistical models
******************

Struct-LMM
^^^^^^^^^^

Struct-LMM can be use to test for interaction with multi-dimensional environments or
to test for association of genetic variants while accounting for GxE interactions.
The Struct-LMM model is

.. math::
    \mathbf{y}=
    \underbrace{\mathbf{F}\mathbf{b}}_{\text{covariates}}+
    \underbrace{\mathbf{g}\beta}_{\text{genetics}}+
    \underbrace{\mathbf{g}\odot\boldsymbol{\gamma}}_{\text{G$\times$E}}+
    \underbrace{\mathbf{u}}_{\text{random effect}}+
    \underbrace{\boldsymbol{\psi}}_{\text{noise}}

where

.. math::
    \boldsymbol{\gamma}\sim\mathcal{N}(\mathbf{0},
    \underbrace{\sigma^2_h\boldsymbol{EE}^T}_{\text{GxE}})

.. math::
    \mathbf{u}\sim\mathcal{N}(\mathbf{0}, \sigma_u^2\mathbf{R}^T)

.. math::
    \boldsymbol{\psi}\sim\mathcal{N}(\mathbf{0}, \sigma_n^2\mathbf{I}_N)

It will be a doctest.

from limix.qtl import GWAS_StructLMM
random = RandomState(1)
envs = random.randn(pheno.shape[0], 30)
slmm = GWAS_StructLMM(pheno, envs, covs=covs, tests=['inter', 'assoc'],
verbose=True)
res = slmm.process(snps[:,:5])
print(res.head())
pvi       pva
0  0.991105  0.926479
1  0.956181  0.984790
2  0.954051  0.989192
3  0.997851  0.393730
4  0.946831  0.375530

The process method returns two sets of P values:
(i) ``pvi`` are the interaction P values,
(ii) ``pva`` are the association P values.
