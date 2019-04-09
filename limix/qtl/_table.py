class Table:
    def __init__(self, columns, index):
        self._columns = list(columns)
        self._index = list(index)
        self._values = None

    def set_values(self, values):
        from numpy import ndarray, atleast_2d

        if isinstance(values, ndarray):
            values = atleast_2d(values.T).T
        if hasattr(values, "tolist"):
            self._values = values.tolist()
        else:
            self._values = list(values)

    def draw(self):
        from texttable import Texttable

        table = Texttable(max_width=88)
        table.set_deco(Texttable.HEADER)
        table.set_chars(["", "", "", "-"])
        table.set_cols_dtype(["t"] + ["e"] * len(self._values[0]))
        table.set_cols_align(["l"] + ["r"] * len(self._values[0]))

        rows = [[""] + self._columns]
        for i, r in enumerate(self._values):
            rows.append([self._index[i]] + list(r))
        table.add_rows(rows)

        return table.draw()
