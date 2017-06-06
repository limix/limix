# TODO for 1.1.0

- [ ] limix.qtl.qtl_test_lmm: 10k by 500k datasets, low memory footprint
      (i.e., <=60 GB).
- [ ] limix.qtl.qtl_test_glmm: 10k by 500k datasets, low memory footprint
      (i.e., <=60 GB).
- [ ] Both limix.qtl.qtl_test_glmm and limix.qtl.qtl_test_lmm should accept
      and take advantage of dask array.
- [ ] limix.qtl.qtl_test_lmm: show progress bar like limix.qtl.qtl_test_glmm.
- [ ] limix.qtl.qtl_test_lmm, limix.qtl.qtl_test_glmm: raise exception for
      wrong parameters (NaN, Inf) and coerce sequences to numpy arrays.
- [ ] Standardize arguments, set labels, set title, set colors, for plotting
      functions: limix.plot.plot_manhattan, limix.plot.qqplot,
      limix.plot.plot_kinship.
- [ ] Define ConvergenceError exception to be raised whenever our methods
      doesn't converge.
- [ ] Convert documentation style to numpy style, of every function.

# TODO for 1.2.0

- [ ] Nicely print a summary of a QTL and heritability result.
- [ ] I/O for BGEN file format.
- [ ] Add full-world QTL and heritability analysis into readthedocs.

# Limix

[![PyPI-License](https://img.shields.io/pypi/l/limix.svg?style=flat-square)](https://pypi.python.org/pypi/limix/)
[![PyPI-Version](https://img.shields.io/pypi/v/limix.svg?style=flat-square)](https://pypi.python.org/pypi/limix/)
[![Documentation Status](https://readthedocs.org/projects/limix/badge/?style=flat-square&version=latest)](https://limix.readthedocs.io/)

Genomic analyses require flexible models that can be adapted to the needs of the user.
Limix is a flexible and efficient linear mixed model library with interfaces to Python.

Limix includes methods for
- single-variant association and interaction testing,
- variance decompostion analysis with linear mixed models,
- association and interaction set tests,
- as well as different utils for statistical analysis, basic i/o and plotting.

A description of the public interface is found at
https://limix.readthedocs.io/.

iPython notebook tutorials are available from github repository:
https://github.com/limix/limix-tutorials.

These tutorials can also be viewed using the ipython notebook viewer:
http://nbviewer.ipython.org/github/limix/limix-tutorials/blob/master/index.ipynb.

#### Highlights

- **iSet**: an interaction set test with external context
  ([paper](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1006693))
  ([public interface](http://limix.readthedocs.io/iSet.html))
  ([tutorial](https://github.com/limix/limix-tutorials/tree/master/iSet))

- **mtSet**: an efficient multi-trait set test
  ([paper](http://www.nature.com/nmeth/journal/v12/n8/abs/nmeth.3439.html))
  ([public interface](http://limix.readthedocs.io/mtSet.html))
  ([tutorial](https://github.com/limix/limix-tutorials/tree/master/mtSet))

## Install

The recommended way of installing it is via
[pip](https://pypi.python.org/pypi/pip)

```bash
pip install limix
```

## Problems

If you encounter any issue, please, [submit it](https://github.com/limix/limix/issues).

## Authors

* **Christoph Lippert** - [https://github.com/clippert](https://github.com/clippert)
* **Danilo Horta** - [https://github.com/Horta](https://github.com/Horta)
* **Francesco Paolo Casale** - [https://github.com/fpcasale](https://github.com/fpcasale)
* **Oliver Stegle** - [https://github.com/ostegle](https://github.com/ostegle)

## License

This project is licensed under the Apache License (Version 2.0, January 2004) -
see the [LICENSE](LICENSE) file for details
