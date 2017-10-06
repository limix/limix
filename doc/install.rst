*******
Install
*******

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

If it does not work, we recommend to consider the first option again.
It only requires the installation of the excelent Python and R platform
`Anaconda`_ designed for reseachers, and then perform a single command.

.. _liknorm: https://github.com/limix/liknorm
.. _conda: http://conda.pydata.org/docs/index.html
.. _pip: https://pypi.python.org/pypi/pip
.. _hcephes: https://github.com/limix/hcephes
.. _Anaconda: https://www.continuum.io/downloads
