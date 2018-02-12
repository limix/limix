import os
import sys

import sphinx_rtd_theme

sys.path.insert(0, os.path.abspath('..'))


def get_version():
    import limix
    return limix.__version__


extensions = [
    'sphinx.ext.autodoc', 'sphinx.ext.doctest', 'sphinx.ext.viewcode',
    'sphinx.ext.intersphinx', 'sphinx.ext.napoleon', 'sphinx.ext.mathjax',
    'matplotlib.sphinxext.only_directives',
    'matplotlib.sphinxext.plot_directive',
    'IPython.sphinxext.ipython_directive',
    'IPython.sphinxext.ipython_console_highlighting'
]

templates_path = ['_templates']

source_suffix = '.rst'

master_doc = 'index'

project = 'limix'
copyright = '2018, Danilo Horta'
author = 'Danilo Horta'

version = get_version()
release = version

language = None

exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', 'conf.py']

pygments_style = 'sphinx'

todo_include_todos = False

html_theme = 'sphinx_rtd_theme'
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

html_sidebars = {
    '**': [
        'relations.html',
        'searchbox.html',
    ]
}

htmlhelp_basename = 'limixdoc'

man_pages = [(master_doc, 'limix', 'limix Documentation', [author], 1)]

epub_title = project
epub_author = author
epub_publisher = author
epub_copyright = copyright

epub_exclude_files = ['search.html']

intersphinx_mapping = {
    'python': ('http://docs.python.org/', None),
    'numpy': ('http://docs.scipy.org/doc/numpy/', None),
}
