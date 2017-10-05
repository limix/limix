
Limix
=====

|PyPI-Status| |Conda-Forge-Status| |Conda-Downloads|

|Build-Status| |Codacy-Grade| |License-Badge| |Doc-Status|


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

Install
-------

The recommended way of installing it is via conda_

.. code:: bash

    conda install -c conda-forge limix

It will handle all the necessary dependencies and should work for GNU/Linux
distributions, MacOS, and Windows.

An alternative way would be via pip_
In this case, first you need to install hcephes_ and liknorm_ C libraries.
Then,

.. code:: bash

    pip install limix


Running the tests
-----------------

After installation, you can test it

.. code:: bash

    python -c "import limix; limix.test()"

as long as you have pytest_.

Authors
-------

* `Christoph Lippert`_
* `Danilo Horta`_
* `Francesco Paolo Casale`_
* `Oliver Stegle`_

License
-------
This project is licensed under the Apache License License - see the
`License file`_ for details.


.. |Build-Status| image:: https://travis-ci.org/limix/limix.svg?branch=refactoring
    :target: https://travis-ci.org/limix/limix

.. |Codacy-Grade| image:: https://api.codacy.com/project/badge/Grade/e0227434c8f040888ff92d1a4d67bcc8
    :target: https://www.codacy.com/app/danilo.horta/limix?utm_source=github.com&utm_medium=referral&utm_content=limix/limix&utm_campaign=badger

.. |PyPI-Status| image:: https://img.shields.io/pypi/v/limix.svg
    :target: https://pypi.python.org/pypi/limix

.. |PyPI-Versions| image:: https://img.shields.io/pypi/pyversions/limix.svg
    :target: https://pypi.python.org/pypi/limix

.. |Conda-Forge-Status| image:: https://anaconda.org/conda-forge/limix/badges/version.svg
    :target: https://anaconda.org/conda-forge/limix

.. |Conda-Downloads| image:: https://anaconda.org/conda-forge/limix/badges/downloads.svg
    :target: https://anaconda.org/conda-forge/limix

.. |License-Badge| image:: https://img.shields.io/pypi/l/limix.svg
    :target: https://raw.githubusercontent.com/limix/limix/refactoring/LICENSE.txt

.. |Doc-Status| image:: https://readthedocs.org/projects/limix/badge/?style=flat-square&version=refactoring
    :target: https://limix.readthedocs.io/

.. |Coverage| image:: https://coveralls.io/repos/github/limix/limix/badge.svg?branch=refactoring
    :target: https://coveralls.io/github/limix/limix?branch=refactoring

.. _License file: https://raw.githubusercontent.com/limix/limix/refactoring/LICENSE.txt

.. _Christoph Lippert: https://github.com/clippert

.. _Danilo Horta: https://github.com/horta

.. _Francesco Paolo Casale: https://github.com/fpcasale

.. _Oliver Stegle: https://github.com/ostegle

.. _conda: http://conda.pydata.org/docs/index.html

.. _pip: https://pypi.python.org/pypi/pip

.. _pytest: http://docs.pytest.org/en/latest/

.. _hcephes: https://github.com/limix/hcephes

.. _liknorm: https://github.com/limix/liknorm
