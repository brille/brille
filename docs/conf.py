# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
import shlex
import subprocess
sys.path.insert(0, os.path.abspath('.'))

# -- Project information -----------------------------------------------------

project = 'brille'
copyright = '2020, Gregory Tucker'
author = 'Gregory Tucker'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.intersphinx',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',
    # 'breathe',
    'sphinxcontrib.katex',
    'sphinxcontrib.tikz',
    # 'exhale',
]

napoleon_use_ivar = True
napoleon_use_param = False
napoleon_use_admonition_for_notes = True

tikz_proc_suite = 'GhostScript'

breathe_projects = {'brille' : '_build/doxygenxml/'}
breathe_default_project = 'brille'
breathe_domain_by_extension = {'h': 'cpp', 'h': 'hpp', 'h': 'tpp'}

# setup the exhale extension
exhale_args = {
    # required arguments first:
    "containmentFolder": "./api",
    "rootFileName": "library_root.rst",
    "rootFileTitle": "C++ Library API",
    "doxygenStripFromPath": "..",
    # Suggested optional arguments
    "createTreeView": True,
    # TIP: if using the sphinx-bootstrap-theme, you also need
    # "treeViewIsBootstrap": True,
    "exhaleExecutesDoxygen": True,
    #"exhaleSilentDoxygen": True,
    "exhaleUseDoxyfile": True,
    # "exhaleDoxygenStdin": "INPUT = ../wrap/", # "INPUT = ../src/",
}

autosummary_generate = True

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

source_suffice = '.rst'
master_doc = 'index'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------
# # following github.com/pybind/pybind11/blob/stable/docs/conf.py:
# this_is_rtd = os.environ.get('READTHEDOCS',None) == 'True'
# if this_is_rtd:
#     html_context = {'css_files' :[
#         '//media.readthedocs.org/css/sphinx_rtd_theme.css',
#         '//media.readthedocs.org/css/readthedocs-doc-embed.css',
#         '_static/theme_overrides.css'
#     ]}
# else:
#     import sphinx_rtd_theme
#     html_theme = 'sphinx_rtd_theme'
#     html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]
#     html_context = {'css_files': ['_static/theme_overrides.css']}

# Alternatively, just specify a theme:
html_theme = 'sphinx_rtd_theme'
html_logo = '../brille.svg'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# again, following from github.com/pybind/pybind11/blob/stable/docs/conf.py:
def generate_doxygen_xml(app):
    build_dir = os.path.join(app.confdir, '_build')
    if not os.path.exists(build_dir):
        os.mkdir(build_dir)

    try:
        subprocess.call(['doxygen', '--version'])
        retcode = subprocess.call(['doxygen'], cwd=app.confdir)
        if retcode < 0:
            sys.stderr.write("doxygen error code: {}\n".format(-retcode))
    except OSError as e:
        sys.stderr.write("doxygen execution failed: {}\n".format(e))

def setup(app):
    """Add hook for building Doxygen xml when needed"""
    #app.connect("builder-inited", generate_doxygen_xml)
