# limix

[![Travis](https://img.shields.io/travis/limix/limix.svg?style=flat-square&label=linux%20%2F%20macos%20build)](https://travis-ci.org/limix/limix) [![AppVeyor](https://img.shields.io/appveyor/ci/Horta/limix.svg?style=flat-square&label=windows%20build)](https://ci.appveyor.com/project/Horta/limix) [![Documentation](https://readthedocs.org/projects/limix/badge/?version=stable&style=flat-square)](https://limix.readthedocs.io/)

Genomic analyses require flexible models that can be adapted to the needs of the user.

Limix is a flexible and efficient linear mixed model library with interfaces to Python.
It includes methods for

- Single-variant association and interaction testing
- Variance decompostion analysis with linear mixed models
- Association and interaction set tests
- Different utils for statistical analysis, basic i/o, and plotting.

The documentation can be found at https://limix.readthedocs.io.

## Install

Installation is easy and works on macOS, Linux, and Windows:

```bash
pip install limix
```

If you already have Limix but want to upgrade it to the latest version:

```bash
pip install limix --upgrade
```

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
