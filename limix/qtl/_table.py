class Table:
    def __init__(self, columns, index):
        self._columns = list(columns)
        self._index = list(index)
        self._values = []

    def add_column(self, values):
        from numpy import atleast_1d

        v = atleast_1d(values).ravel().tolist()
        self._values.append(v)

    def draw(self):
        from texttable import Texttable

        table = Texttable(max_width=88)
        table.set_deco(Texttable.HEADER)
        table.set_chars(["", "", "", "-"])
        ncols = len(self._values)
        table.set_cols_dtype(["t"] + ["e"] * ncols)
        table.set_cols_align(["l"] + ["r"] * ncols)

        rows = [[""] + self._columns]
        for i, r in enumerate(zip(*self._values)):
            rows.append([self._index[i]] + list(r))

        table.add_rows(rows)

        return table.draw()
