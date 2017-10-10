def call(cmd):
    r"""Execute a shell command.

    Parameters
    ----------
    cmd : str
        Command to be executed.

    Returns
    -------
    int
        Exit status.
    """
    import subprocess

    return subprocess.call(cmd, shell=True)
