from setuptools import setup

if __name__ == "__main__":
    console_scripts = ["limix = limix.cmdlimix:entry"]
    setup(
        dependency_links=[
            "git+https://github.com/horta/pytest-doctestplus.git@precision#egg=pytest-doctestplus-0.1.4.dev"
        ],
        entry_points=dict(console_scripts=console_scripts),
    )
