class Aligned:
    def __init__(self):
        self._items = []

    def add_item(self, field, value):
        self._items.append([field, value])

    def draw(self):
        from textwrap import TextWrapper

        max_len = max(len(i[0]) for i in self._items)
        max_len += 3

        for i in self._items:
            i[0] = i[0] + " " * (max_len - 2 - len(i[0]))
            i[0] += "= "

        s = " " * max_len

        msg = ""
        for i in self._items:
            wrapper = TextWrapper(initial_indent=i[0], width=88, subsequent_indent=s)
            msg += wrapper.fill(str(i[1])) + "\n"

        return msg
