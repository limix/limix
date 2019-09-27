def draw_title(title):
    msg = f"{title}\n"
    msg += "-" * len(title) + "\n"
    return msg


def draw_dataframe(title, df):
    msg = repr(df)
    k = msg.find("\n") - len(title) - 2
    left = ("-" * (k // 2)) + " "
    right = " " + ("-" * (k // 2 + k % 2))
    out = left + title + right + "\n"
    out += msg
    return out


def draw_list(x, n):
    x = list(x)
    if len(x) <= n:
        return str(x)
    if n < 1:
        raise ValueError("`n` must be greater than `0` for non-empty lists.")
    if n == 1:
        return "[...]"
    if n == 2:
        return "[..., " + str(x[-1:])[1:]
    return str(x[: n - 2])[:-1] + ", ..., " + str(x[-1:])[1:]
