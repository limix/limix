***********************
Heritability estimation
***********************

We provide heritability estimation for Normal, Bernoulli, Probit, Binomial, and Poisson
phenotypes.
A standard LMM is used for Normal traits:

.. math::

    𝐲 = 𝙼𝛂 + 𝐯 + 𝛆,

where

.. math::

    𝐯 ∼ 𝓝(𝟎, 𝓋₀𝙺) ~~\text{and}~~ 𝛆 ∼ 𝓝(𝟎, 𝓋₁𝙸).

A GLMM is used to model the other type of traits:

.. math::

    𝐳 = 𝙼𝛂 + 𝐯 + 𝛆, ~~\text{where}~~ yᵢ|𝐳 ∼ 𝙴𝚡𝚙𝙵𝚊𝚖(𝜇ᵢ=g(zᵢ))

and 𝐯 and 𝛆 are defined as before.

In both cases, the parameters are the same: 𝛂, 𝓋₀, and 𝓋₁. They are fitted via
restricted maximum likelihood for LMM and via maximum likelihood for GLMM.
The covariance-matrix 𝙺 given by the user is normalised before the model is fitted as
follows:

.. code-block:: python

    K = K / K.diagonal().mean()

.. autofunction:: limix.her.estimate
    :noindex:
