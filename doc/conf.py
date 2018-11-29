import errno
import os
import shutil
import sys
import time
from os.path import join

import sphinx_rtd_theme

sys.path.insert(0, os.path.abspath(".."))


def _get_version():
    import limix

    return limix.__version__


def _get_name():
    import limix

    return limix.__name__


extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.doctest",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx.ext.napoleon",
    "sphinx.ext.mathjax",
    "sphinxcontrib.programoutput",
    "matplotlib.sphinxext.plot_directive",
    "IPython.sphinxext.ipython_directive",
    "IPython.sphinxext.ipython_console_highlighting",
    "sphinx.ext.autosummary",
]

autodoc_default_flags = ["members"]
autodoc_mock_imports = ["_tkinter"]
autosummary_generate = True
napoleon_numpy_docstring = True
templates_path = ["_templates"]

source_suffix = ".rst"

master_doc = "index"

project = _get_name()
copyright = "2018, Danilo Horta"
author = "Danilo Horta"

version = _get_version()
release = version
today = time.strftime("%B %d, %Y")

language = None

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "conf.py"]

pygments_style = "default"

todo_include_todos = False

html_theme = "sphinx_rtd_theme"
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

html_sidebars = {"**": ["relations.html", "searchbox.html"]}
htmlhelp_basename = "{}doc".format(project)

man_pages = [(master_doc, _get_name(), "{} documentation".format(project), [author], 1)]

intersphinx_mapping = {
    "python": ("http://docs.python.org/", None),
    "numpy": ("https://docs.scipy.org/doc/numpy/", None),
    "limix-plot": ("https://limix-plot.readthedocs.io/en/stable/", None),
    "pandas": ("https://pandas.pydata.org/pandas-docs/stable/", None),
    "dask": ("http://dask.pydata.org/en/latest/", None),
    "xarray": ("http://xarray.pydata.org/en/latest/", None),
}

plot_working_directory = "_build"


def makedirs(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


dir_path = os.path.dirname(os.path.realpath(__file__))

makedirs(join(dir_path, "_build/"))
shutil.copy2(join(dir_path, "matplotlibrc"), join(dir_path, "_build/"))

plot_formats = [("png", 80)]
