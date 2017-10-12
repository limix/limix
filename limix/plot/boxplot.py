def boxplot(df, style=None):
    r"""Box-plot a data frame.

    Parameters
    ----------
    df : data_frame
        A data frame containing `value`, `variable, and `category` columns.
    style : dict
        Keyword arguments forwarded to :func:`matplotlib.axes.Axes.plot`
        function.

    Examples
    --------
    .. plot::

        import seaborn as sns
        import limix

        df = sns.load_dataset("exercise")
        df.rename(
            columns=dict(time='category', kind='variable', pulse='value'),
            inplace=True)

        limix.plot.set_paper_style()
        limix.plot.boxplot(df)
        limix.plot.show()
    """
    from seaborn import boxplot

    if style is None:
        style = dict()

    if 'flierprops' not in style:
        style['flierprops'] = dict()

    style['flierprops'].update(dict(markersize=3.0))

    g = boxplot(y='value', x='category', hue='variable', data=df, **style)

    g.axes.grid(
        True,
        which='major',
        axis='y',
        linewidth=0.75,
        linestyle='-',
        color='#EEEEEE',
        alpha=1.0)

    g.axes.legend().remove()

    return g
