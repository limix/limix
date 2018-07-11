from __future__ import unicode_literals

import os
import sys

from setuptools import find_packages, setup

try:
    import pypandoc

    long_description = pypandoc.convert_file("README.md", "rst")
except (OSError, IOError, ImportError):
    long_description = open("README.md").read()


def setup_package():
    src_path = os.path.dirname(os.path.abspath(sys.argv[0]))
    old_path = os.getcwd()
    os.chdir(src_path)
    sys.path.insert(0, src_path)

    needs_pytest = {"pytest", "test", "ptr"}.intersection(sys.argv)
    pytest_runner = ["pytest-runner>=2.9"] if needs_pytest else []

    setup_requires = ["cython", "numpy>=1.14.5"] + pytest_runner
    install_requires = [
        "asciitree>=0.3.3",
        "dask[array,bag,dataframe,delayed]>=0.18.1",
        "distributed",
        "glimix-core>=1.5.0",
        "h5py",
        "joblib>=0.11",
        "limix-core>=1.0.1",
        "limix-legacy>=0.8.12",
        "ncephes>=1.0.40",
        "numpy-sugar>=1.2.5",
        "numpy>=1.14.5",
        "pandas-plink>=1.2.26",
        "scikit-learn",
        "scipy>=1.1.0",
        "tqdm>=4.10",
    ]
    tests_require = ["pytest", "pytest-console-scripts", "pytest-pep8"]

    console_scripts = [
        "limix_runner=limix.scripts.limix_runner:entry_point",
        "mtset_postprocess=limix.scripts.mtset_postprocess:entry_point",
        "mtset_preprocess=limix.scripts.mtset_preprocess:entry_point",
        "mtset_definesets=limix.scripts.mtset_definesets:entry_point",
        "mtset_simpheno=limix.scripts.mtset_simpheno:entry_point",
        "mtset_analyze=limix.scripts.mtset_analyze:entry_point",
        "limix_converter=limix.scripts.limix_converter:entry_point",
        "iset_analyze=limix.scripts.iset_analyze:entry_point",
        "iset_postprocess=limix.scripts.iset_postprocess:entry_point",
        "limix=limix.scripts.limix:entry_point",
        "ilimix=limix.scripts.ilimix:entry_point",
    ]

    metadata = dict(
        name="limix",
        version="1.0.18",
        maintainer="Limix Developers",
        maintainer_email="horta@ebi.ac.uk",
        author=(
            "Christoph Lippert, Danilo Horta, "
            + "Francesco Paolo Casale, Oliver Stegle"
        ),
        author_email="stegle@ebi.ac.uk",
        license="Apache License 2.0",
        description="A flexible and fast mixed model toolbox.",
        long_description=long_description,
        url="https://github.com/limix/limix",
        packages=find_packages(),
        zip_safe=False,
        install_requires=install_requires,
        setup_requires=setup_requires,
        tests_require=tests_require,
        include_package_data=True,
        classifiers=[
            "Development Status :: 5 - Production/Stable",
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent",
        ],
        entry_points={"console_scripts": console_scripts},
    )

    try:
        setup(**metadata)
    finally:
        del sys.path[0]
        os.chdir(old_path)


if __name__ == "__main__":
    setup_package()
