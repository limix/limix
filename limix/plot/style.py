def set_paper_style():
    from seaborn import set_context, despine, set_style

    set_context(
        "paper", rc={"font.size": 8,
                     "axes.titlesize": 8,
                     "axes.labelsize": 5})
    despine(left=False, right=False, top=False, bottom=False)

    set_style({
        'font.sans-serif': [
            'Helvetica', 'Arial', 'Baskerville', 'Caslon', 'Garamond',
            'DejaVu Sans', 'Bitstream Vera Sans', 'Computer Modern Sans Serif',
            'Lucida Grande', 'Verdana', 'Geneva', 'Lucid', 'Avant Garde',
            'sans-serif'
        ]
    })
