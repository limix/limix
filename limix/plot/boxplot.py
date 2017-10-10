def boxplot(df, style=None):
    from seaborn import factorplot

    g = factorplot(
        kind='box',
        y='value',
        x='category',
        hue='variable',
        data=df,
        legend=False)

    g.axes[0, 0].grid(True, which='major', axis='y', color="#DDDDDD")

    return g


if __name__ == '__main__':
    import seaborn as sns
    import limix

    df = sns.load_dataset("exercise")
    df.rename(
        columns=dict(time='category', kind='variable', pulse='value'),
        inplace=True)

    limix.plot.boxplot(df)
    limix.plot.set_paper_style()

    limix.plot.show()
