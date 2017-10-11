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

        limix.plot.boxplot(df)
        limix.plot.set_paper_style()
        limix.plot.show()
    """
    from seaborn import factorplot

    if style is None:
        style = {label: dict() for label in labels}

    g = factorplot(
        kind='box',
        y='value',
        x='category',
        hue='variable',
        data=df,
        legend=False,
        **style)

    g.axes[0, 0].grid(
        True, which='major', axis='y', color="#AAAAAA", alpha=0.2)

    return g
