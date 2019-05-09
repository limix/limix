# limix

[![Travis](https://img.shields.io/travis/com/limix/limix.svg?style=flat-square&label=linux%20%2F%20macos%20build)](https://travis-ci.com/limix/limix) [![AppVeyor](https://img.shields.io/appveyor/ci/Horta/limix.svg?style=flat-square&label=windows%20build)](https://ci.appveyor.com/project/Horta/limix) [![Documentation](https://readthedocs.org/projects/limix/badge/?version=latest&style=flat-square)](https://limix.readthedocs.io/)

Genomic analyses require flexible models that can be adapted to the needs of the user.

Limix is a flexible and efficient linear mixed model library with interfaces to Python.
It includes methods for

- Single-variant association and interaction testing
- Variance decompostion analysis with linear mixed models
- Association and interaction set tests
- Different utils for statistical analysis, basic i/o, and plotting.

The documentation can be found at https://limix.readthedocs.io.

## Install

> **NOTE**: We will be maintaining limix 2.0.x for a while, in case you find some
> missing feature in limix 3.0.x. If that is the case, please, type `pip install "limix <3,>=2"` in your terminal.

Installation is easy and works on macOS, Linux, and Windows:

```bash
pip install limix
```

If you already have Limix but want to upgrade it to the latest version:

```bash
pip install limix --upgrade
```

## Interactive tutorials

- [eQTL](https://mybinder.org/v2/gh/limix/limix-tutorials/master?filepath=eQTL.ipynb) (requires limix 2.0.x)
- [Struct-LMM](https://mybinder.org/v2/gh/limix/limix-tutorials/master?filepath=struct-lmm.ipynb) (requires limix 2.0.x)

## Running tests

After installation, you can test it

```bash
python -c "import limix; limix.test()"
```

as long as you have [pytest](https://docs.pytest.org/en/latest/).

## Authors

* [Christoph Lippert](https://github.com/clippert)
* [Danilo Horta](https://github.com/horta)
* [Francesco Paolo Casale](https://github.com/fpcasale)
* [Oliver Stegle](https://github.com/ostegle)

## License

This project is licensed under the [Apache License License](https://raw.githubusercontent.com/limix/limix/2.0.0/LICENSE.md).
