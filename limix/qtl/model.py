from __future__ import division

from pandas import DataFrame, Series

from limix.stats import effsizes_se, lrt_pvalues


class QTLModel(object):
    r"""Result of a QTL analysis.

    An instance of this class is returned by :func:`limix.qtl.scan`.
    """

    def __init__(self, null_lml, alt_lmls, effsizes, null_covariate_effsizes):
        self._null_lml = null_lml
        self._alt_lmls = alt_lmls
        self._effsizes = effsizes
        self._null_covariate_effsizes = null_covariate_effsizes

    @property
    def null_lml(self):
        r"""Log of the marginal likelihood under the null hypothesis.

        Returns
        -------
        float
            Log of marginal likelihood.
        """
        return self._null_lml

    @property
    def alt_lmls(self):
        r"""Log of the marginal likelihoods across tested variants.

        Returns
        -------
        array_like
            Log of marginal likelihoods.
        """
        return self._alt_lmls

    @property
    def variant_effsizes(self):
        r"""Variant effect-sizes.

        Returns
        -------
        array_like
            Estimated variant effect sizes.
        """
        return self._effsizes

    @property
    def variant_pvalues(self):
        r"""Variant p-values.

        Returns
        -------
        array_like
            Association significance between variant and phenotype.
        """
        pv = lrt_pvalues(self.null_lml, self.alt_lmls.values)
        return Series(pv, list(self.alt_lmls.keys()))

    @property
    def variant_effsizes_se(self):
        r"""Standard errors of the variant effect sizes.

        Returns
        -------
        array_like
            Estimated standard errors of the variant effect sizes.
        """
        ese = effsizes_se(self.variant_effsizes.values,
                          self.variant_pvalues.values)
        return Series(ese, list(self.alt_lmls.keys()))

    @property
    def null_covariate_effsizes(self):
        r"""Covariate effect-sizes under the null hypothesis.

        Returns
        -------
        array_like
            Estimated covariant effect sizes under the null hypothesis.
        """
        return self._null_covariate_effsizes

    def __str__(self):
        data = dict(
            effsizes=self.variant_effsizes,
            effsizes_se=self.variant_effsizes_se,
            pvalues=self.variant_pvalues)

        variant_msg = str(DataFrame(data=data).describe())

        data = self.null_covariate_effsizes
        k = data.index.values
        v = [[vi] for vi in data.values]

        covariate_msg = str(DataFrame(data=dict(zip(list(k), list(v)))))
        covariate_msg = '\n'.join([x[2:] for x in covariate_msg.split('\n')])

        msg = 'Variants\n' + variant_msg
        msg += '\n\nCovariate effect sizes for the'
        msg += ' null model\n' + covariate_msg

        return msg


class IQTLModel(object):
    r"""Result of a Interaction QTL analysis.

    An instance of this class is returned by :func:`limix.qtl.iscan`.
    """

    def __init__(self, null_lml, alt_lmls, effsizes, null_covariate_effsizes):
        self._null_lml = null_lml
        self._alt_lmls = alt_lmls
        self._effsizes = effsizes
        self._null_covariate_effsizes = null_covariate_effsizes

    @property
    def null_lml(self):
        r"""Log of the marginal likelihood under the null hypothesis.

        Returns
        -------
        float
            Log of marginal likelihood.
        """
        return self._null_lml

    @property
    def alt_lmls(self):
        r"""Log of the marginal likelihoods across tested variants.

        Returns
        -------
        array_like
            Log of marginal likelihoods.
        """
        return self._alt_lmls

    @property
    def variant_effsizes(self):
        r"""Variant effect-sizes.

        Returns
        -------
        array_like
            Estimated variant effect sizes.
        """
        return self._effsizes

    @property
    def variant_pvalues(self):
        r"""Variant p-values.

        Returns
        -------
        array_like
            Association significance between variant and phenotype.
        """
        pv = dict()

        for (i, null_lml) in enumerate(self.null_lml):
            alt_lmls = self.alt_lmls.iloc[i, :].values
            pv[self.null_lml.index[i]] = lrt_pvalues(null_lml, alt_lmls)

        pv = DataFrame(data=pv, index=list(self.alt_lmls.columns))
        return pv.transpose()

    @property
    def variant_effsizes_se(self):
        r"""Standard errors of the variant effect sizes.

        Returns
        -------
        array_like
            Estimated standard errors of the variant effect sizes.
        """
        ese = dict()

        for i, ve in self.variant_effsizes.iterrows():
            vpvalues = self.variant_pvalues.loc[i, :].values
            ese[i] = effsizes_se(ve.values, vpvalues)

        ese = DataFrame(data=ese, index=list(self.variant_pvalues.columns))
        return ese.transpose()

    @property
    def null_covariate_effsizes(self):
        r"""Covariate effect-sizes under the null hypothesis.

        Returns
        -------
        array_like
            Estimated covariant effect sizes under the null hypothesis.
        """
        return self._null_covariate_effsizes

    # def __str__(self):
    #     from pandas import DataFrame
    #
    #     data = dict(
    #         effsizes=self.variant_effsizes,
    #         effsizes_se=self.variant_effsizes_se,
    #         pvalues=self.variant_pvalues)
    #
    #     variant_msg = str(DataFrame(data=data).describe())
    #
    #     data = self.null_covariate_effsizes
    #     k = data.index.values
    #     v = [[vi] for vi in data.values]
    #
    #     covariate_msg = str(DataFrame(data=dict(zip(list(k), list(v)))))
    #     covariate_msg = '\n'.join([x[2:] for x in covariate_msg.split('\n')])
    #
    #     msg = 'Variants\n' + variant_msg
    #     msg += '\n\nCovariate effect sizes for the'
    #     msg += ' null model\n' + covariate_msg
    #
    #     return msg
