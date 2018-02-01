# Limix

Genomic analyses require flexible models that can be adapted to the needs of the user.
Limix is a flexible and efficient linear mixed model library with interfaces to Python.

Limix includes methods for

- Single-variant association and interaction testing
- Variance decompostion analysis with linear mixed models
- Association and interaction set tests
- Different utils for statistical analysis, basic i/o, and plotting.

The documentation can be found at  https://limix.readthedocs.io/en/develop.

Install
-------

The development version of limix can be installed on MacOS and Linux via

.. code:: bash

    bash <(curl -fsSL https://raw.githubusercontent.com/limix/limix/develop/install)

Stable versions of limix are installed via conda_ though

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


.. |Build-Status| image:: https://travis-ci.org/limix/limix.svg?branch=develop
    :target: https://travis-ci.org/limix/limix

.. |Codacy-Grade| image:: https://api.codacy.com/project/badge/Grade/cd0ff739fa004091a1459f1a13d55ad0
    :target: https://www.codacy.com/app/danilo.horta/limix?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=limix/limix&amp;utm_campaign=Badge_Grade

.. |PyPI-Status| image:: https://img.shields.io/pypi/v/limix.svg
    :target: https://pypi.python.org/pypi/limix

.. |PyPI-Versions| image:: https://img.shields.io/pypi/pyversions/limix.svg
    :target: https://pypi.python.org/pypi/limix

.. |Conda-Forge-Status| image:: https://anaconda.org/conda-forge/limix/badges/version.svg
    :target: https://anaconda.org/conda-forge/limix

.. |Conda-Downloads| image:: https://anaconda.org/conda-forge/limix/badges/downloads.svg
    :target: https://anaconda.org/conda-forge/limix

.. |License-Badge| image:: https://img.shields.io/pypi/l/limix.svg
    :target: https://raw.githubusercontent.com/limix/limix/develop/LICENSE.txt

.. |Doc-Status| image:: https://readthedocs.org/projects/limix/badge/?style=flat-square&version=develop
    :target: https://limix.readthedocs.io/en/develop

.. |Coverage| image:: https://codecov.io/gh/limix/limix/branch/develop/graph/badge.svg
    :target: https://codecov.io/gh/limix/limix/branch/develop

.. _License file: https://raw.githubusercontent.com/limix/limix/develop/LICENSE.txt

.. _Christoph Lippert: https://github.com/clippert

.. _Danilo Horta: https://github.com/horta

.. _Francesco Paolo Casale: https://github.com/fpcasale

.. _Oliver Stegle: https://github.com/ostegle

.. _conda: http://conda.pydata.org/docs/index.html

.. _pip: https://pypi.python.org/pypi/pip

.. _pytest: http://docs.pytest.org/en/latest/

.. _hcephes: https://github.com/limix/hcephes

.. _liknorm: https://github.com/limix/liknorm
