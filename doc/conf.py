import os
import sys

import sphinx_rtd_theme

sys.path.insert(0, os.path.abspath('..'))


def _get_version():
    import limix
    return limix.__version__


def _get_name():
    import glimix_core
    return glimix_core.__name__


extensions = [
    'sphinx.ext.autodoc', 'sphinx.ext.doctest', 'sphinx.ext.viewcode',
    'sphinx.ext.intersphinx', 'sphinx.ext.napoleon', 'sphinx.ext.mathjax',
    'matplotlib.sphinxext.only_directives',
    'matplotlib.sphinxext.plot_directive',
    'IPython.sphinxext.ipython_directive',
    'IPython.sphinxext.ipython_console_highlighting', 'nb2plots',
    'sphinx.ext.autosummary'
]

autodoc_default_flags = ['members']
autosummary_generate = True
napoleon_numpy_docstring = True
templates_path = ['_templates']

source_suffix = '.rst'

master_doc = 'index'

project = _get_name()
copyright = '2018, Danilo Horta'
author = 'Danilo Horta'

version = _get_version()
release = version

language = None

exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', 'conf.py']

pygments_style = 'default'

todo_include_todos = False

html_theme = 'sphinx_rtd_theme'
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

html_sidebars = {
    '**': [
        'relations.html',
        'searchbox.html',
    ]
}
htmlhelp_basename = '{}doc'.format(project)

man_pages = [(master_doc, _get_name(), '{} documentation'.format(project),
              [author], 1)]

intersphinx_mapping = {
    'python': ('http://docs.python.org/', None),
    'limix-plot': ('https://limix-plot.readthedocs.io/en/stable/', None)
}
