*******
Install
*******

The recommended way of installing it is via `conda`_:

.. code-block:: bash

    conda install -c conda-forge limix

It will handle all the necessary dependencies and should work for GNU/Linux
distributions, macOS, and 64-bit Windows.

An alternative way would be via `pip`_.
First, make sure you have at least the following Python packages:

- `NumPy`_
- `SciPy`_
- `Matplotlib`_
- `Setuptools`_
- `Scikit-Learn`_

You also need to install an C library called `liknorm`_, and have a
C/C++ compiler.
Finally, run

.. code-block:: bash

    pip install limix

to perform the installation.
If it does not work, we recommend to consider the first option again.
It only requires the installation of the excelent Python and R platform
`Anaconda`_ designed for reseachers, and then perform a single command.

.. _liknorm: https://github.com/limix/liknorm
.. _conda: http://conda.pydata.org/docs/index.html
.. _pip: https://pypi.python.org/pypi/pip
.. _NumPy: http://www.numpy.org
.. _SciPy: https://www.scipy.org
.. _Matplotlib: https://matplotlib.org
.. _Setuptools: https://pypi.python.org/pypi/setuptools
.. _Scikit-Learn: http://scikit-learn.org
.. _Anaconda: https://www.continuum.io/downloads
