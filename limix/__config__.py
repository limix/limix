def get_info(name):
    import os

    info = {
        "interactive": _interactive_session(),
        "rich_text": _in_terminal() or _in_ipython_frontend(),
        "building_doc": os.environ.get("LIMIX_DOC", "False") == "True",
    }
    return info[name]


def _interactive_session():
    """
    Check if we're running in an interactive shell.

    Returns
    -------
    ``True`` if running under python/ipython interactive shell; ``False`` otherwise.
    """

    def check_main():
        import __main__ as main

        return not hasattr(main, "__file__")

    try:
        return hasattr(globals()["__builtins__"], "__IPYTHON__") or check_main()
    except Exception:
        return check_main()


def _in_ipython_frontend():
    """
    Check if we're inside an an IPython zmq frontend.
    """
    if "get_ipython" not in globals():
        return False
    ip = globals()["get_ipython"]()
    return "zmq" in str(type(ip)).lower()


def _in_terminal():
    """
    Detect if Python is running in a terminal.

    Returns
    -------
    bool
        ``True`` if Python is running in a terminal; ``False`` otherwise.
    """
    # Assume standard Python interpreter in a terminal.
    if "get_ipython" not in globals():
        return True
    ip = globals()["get_ipython"]()
    # IPython as a Jupyter kernel.
    if hasattr(ip, "kernel"):
        return False
    return True
