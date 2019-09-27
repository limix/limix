class AlignedText:
    def __init__(self, sep="= "):
        self._items = []
        self._sep = sep

    def add_item(self, field, value):
        self._items.append([field, value])

    def draw(self):
        from limix import config

        from textwrap import TextWrapper

        max_len = max(len(i[0]) for i in self._items)
        max_len += len(self._sep) + 1

        for i in self._items:
            i[0] = i[0] + " " * (max_len - len(self._sep) - len(i[0]))
            i[0] += self._sep

        s = " " * max_len

        msg = ""
        width = config["display.text_width"]
        for i in self._items:
            wrapper = TextWrapper(initial_indent=i[0], width=width, subsequent_indent=s)
            msg += wrapper.fill(str(i[1])) + "\n"

        return msg.rstrip()
