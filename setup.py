from setuptools import setup

if __name__ == "__main__":
    console_scripts = ["limix = limix.cmdlimix:entry"]
    setup(
        dependency_links=[
            "http://github.com/horta/pytest-doctestplus/tarball/precision#egg=pytest-doctestplus-0.1.4.dev0"
        ],
        entry_points=dict(console_scripts=console_scripts),
    )
