def plot_boxplot(df, style=None, ax=None):
    from seaborn import boxplot

    if style is None:
        style = dict()

    boxplot(
        y='value',
        x='category',
        hue='variable',
        data=df,
        flierprops=dict(markersize=3.0),
        **style)

    ax.grid(
        True,
        which='major',
        axis='y',
        linewidth=0.75,
        linestyle='-',
        color='#EEEEEE',
        alpha=1.0)
