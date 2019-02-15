class Pipeline(object):
    def __init__(self, data):
        self._process = []
        self._data = data
        self._layout = _LayoutChange()

    def append(self, process, name, *args, **kwargs):
        self._process.append(
            {"func": process, "name": name, "args": args, "kwargs": kwargs}
        )

    def run(self, verbose=True):
        from limix._data import conform_dataset, CONF

        target = CONF["varname_to_target"]
        self._layout.append(
            "initial", {target[vn]: self._data[vn].shape for vn in self._data.keys()}
        )

        for p in self._process:
            self._data = p["func"](self._data, *p["args"], **p["kwargs"])
            self._data = conform_dataset(**self._data)

            self._layout.append(
                p["name"],
                {
                    target[n]: self._data[n].shape
                    for n, v in self._data.items()
                    if v is not None
                },
            )

            if self._get_samples().size == 0:
                print(self._layout.to_string())
                raise RuntimeError("Exiting early because there is no sample left.")

        if verbose:
            print(self._layout.to_string())

        return self._data

    def _get_samples(self):
        for x in self._data.values():
            if hasattr(x, "sample"):
                return x.sample
        raise RuntimeError("Could not get samples.")


class _LayoutChange(object):
    def __init__(self):
        self._steps = []

    def append(self, name, target_shapes):
        self._steps.append((name, target_shapes))

    def to_string(self):
        from texttable import Texttable

        table = Texttable()
        header = [""]
        shapes = {k: [k] for k in self._steps[0][1].keys()}

        for step in self._steps:
            header.append(step[0])
            for t, shape in step[1].items():
                shapes[t].append(shape)

        table.header(header)

        table.set_cols_dtype(["t"] * len(header))
        table.set_cols_align(["l"] * len(header))
        table.set_deco(Texttable.HEADER)

        for target in self._steps[0][1].keys():
            table.add_row(shapes[target])

        msg = table.draw()

        msg = self._add_caption(msg, "-", "Table: Data layout transformation.")
        return msg

    def _add_caption(self, msg, c, caption):
        n = len(msg.split("\n")[-1])
        msg += "\n" + (c * n)
        msg += "\n" + caption
        return msg
