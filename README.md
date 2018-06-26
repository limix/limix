# limix

[![Travis](https://img.shields.io/travis/limix/limix.svg?style=flat-square&label=linux%20%2F%20macos%20build)](https://travis-ci.org/limix/limix) [![AppVeyor](https://img.shields.io/appveyor/ci/Horta/limix.svg?style=flat-square&label=windows%20build)](https://ci.appveyor.com/project/Horta/limix) [![Documentation](https://img.shields.io/readthedocs/limix.svg?style=flat-square&version=stable)](https://limix.readthedocs.io/)

Genomic analyses require flexible models that can be adapted to the needs of the user.

Limix is a flexible and efficient linear mixed model library with interfaces to Python.
It includes methods for

- Single-variant association and interaction testing
- Variance decompostion analysis with linear mixed models
- Association and interaction set tests
- Different utils for statistical analysis, basic i/o, and plotting.

The documentation can be found at  https://limix.readthedocs.io/en/develop.

## Install

The development version of limix can be installed on macOS and Linux via

```bash
bash <(curl -fsSL https://raw.githubusercontent.com/limix/limix/2.0.0/install)
```

Stable versions of limix are installed via [conda](http://conda.pydata.org/docs/index.html) though

```bash
conda install -c conda-forge limix
```

It will handle all the necessary dependencies and should work for GNU/Linux
distributions, MacOS, and Windows.

An alternative way would be via [pip](https://pypi.python.org/pypi/pip)
In this case, first you need to install [hcephes](https://github.com/limix/hcephes) and [liknorm](https://github.com/limix/liknorm) C libraries.
Then,

```bash
pip install limix
```


## Running the tests

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
