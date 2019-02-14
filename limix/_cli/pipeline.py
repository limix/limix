from limix._data import CONF


class Pipeline(object):
    def __init__(self, data):
        self._process = []
        self._data = data
        self._layout = _LayoutChange()

    def append(self, process, *args, **kwargs):
        self._process.append({"func": process, "args": args, "kwargs": kwargs})

    def run(self, verbose=True):
        for varname, d in self._data.items():
            target = CONF["varname_to_target"][varname]
            self._layout.append(target, "initial", d.shape)

        for p in self._process:
            self._data = p["func"](self._data, self._layout, *p["args"], **p["kwargs"])

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
        self._targets = {}
        self._steps = ["sentinel"]

    def append(self, target, step, shape):
        if target not in CONF["targets"]:
            breakpoint()
            raise ValueError("Invalid target `{}`.".format(target))

        if target not in self._targets:
            self._targets[target] = {}

        self._targets[target][step] = shape
        if step != self._steps[-1]:
            self._steps.append(step)

    def to_string(self):
        from texttable import Texttable

        table = Texttable()
        header = [""]
        shapes = {k: [k] for k in self._targets.keys()}

        for step in self._steps[1:]:
            header.append(step)
            for target in self._targets.keys():
                v = str(self._targets[target].get(step, "..."))
                shapes[target].append(v)

        table.header(header)

        table.set_cols_dtype(["t"] * len(header))
        table.set_cols_align(["l"] * len(header))
        table.set_deco(Texttable.HEADER)

        for target in self._targets.keys():
            table.add_row(shapes[target])

        msg = table.draw()

        msg = self._add_caption(msg, "-", "Table: Data layout transformation.")
        return msg

    def _add_caption(self, msg, c, caption):
        n = len(msg.split("\n")[-1])
        msg += "\n" + (c * n)
        msg += "\n" + caption
        return msg
