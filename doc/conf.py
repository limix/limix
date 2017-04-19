from __future__ import unicode_literals

import os

try:
    import limix
    version = limix.__version__
except ImportError:
    version = 'unknown'

os.environ['BOKEH_DOCS_MISSING_API_KEY_OK'] = "yes"

extensions = [
    'matplotlib.sphinxext.only_directives',
    'matplotlib.sphinxext.plot_directive',
    'sphinx.ext.autodoc',
    'sphinx.ext.doctest',
    'sphinx.ext.intersphinx',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
]
napoleon_google_docstring = True
templates_path = ['_templates']
source_suffix = '.rst'
master_doc = 'index'
project = 'limix'
copyright = '2017, Christoph Lippert, Danilo Horta, Francesco Paolo Casale, Oliver Stegle'
author = 'Christoph Lippert, Danilo Horta, Francesco Paolo Casale, Oliver Stegle'
release = version
language = None
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', 'conf.py']
pygments_style = 'sphinx'
todo_include_todos = False
html_theme = 'default'
htmlhelp_basename = 'limixdoc'
latex_elements = {}
latex_documents = [
    (master_doc, 'limix.tex', 'limix Documentation',
     'Christoph Lippert, Danilo Horta, Francesco Paolo Casale, Oliver Stegle',
     'manual'),
]
man_pages = [(master_doc, 'limix', 'limix Documentation', [author], 1)]
texinfo_documents = [
    (master_doc, 'limix', 'limix Documentation', author, 'limix',
     'A flexible and fast mixed model toolbox.', 'Miscellaneous'),
]
intersphinx_mapping = {'https://docs.python.org/': None,
                       'http://matplotlib.org': None}
