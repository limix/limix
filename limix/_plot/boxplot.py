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

    ax.grid(True, which='major', axis='y', alpha=1.0)

    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
