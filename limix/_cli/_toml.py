from pathlib import Path


def write_limix_toml(filepath: Path):
    from limix import __version__
    import sys
    from os import getcwd
    import qtoml

    major = sys.version_info.major
    minor = sys.version_info.minor
    micro = sys.version_info.micro
    pyver = f"{major}.{minor}.{micro}"

    with open(filepath / "limix.toml", "w") as f:
        f.write(f"# Limix  {__version__}\n")
        f.write(f"# Python {pyver}\n\n")
        qtoml.dump({"workdir": getcwd(), "cmdline": sys.argv}, f)
