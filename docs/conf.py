"""
Configuration file for the Sphinx documentation builder.

This file only contains a selection of the most common options. For a full
list see the documentation:
https://www.sphinx-doc.org/en/master/usage/configuration.html
"""

# Standard library imports
import os
import sys
from datetime import date

# Sphinx imports
from sphinx_pyproject import SphinxConfig

# Local imports
import kkcalc2

path = os.path.abspath(os.path.dirname(__file__))
pyproj = os.path.join(path, "../pyproject.toml")

config = SphinxConfig(pyproj, globalns=globals())
# Now, variables like 'project', 'version', 'author' will be available
# from the loaded pyproject.toml data.
# You can also access other values through the 'config' object.
# For example: html_theme = config.html_theme if set in pyproject.toml

# -- Options for HTML output ----------------------------------------------
html_theme = config["html_theme"]
html_theme_options = (
    config["html_theme_options"] if "html_theme_options" in config else {}
)
html_sidebars = config["html_sidebars"] if "html_sidebars" in config else {"**": []}
html_context = config["html_context"] if "html_context" in config else {}
html_last_updated_fmt = "%b %d, %Y"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["source/_static"]  # ['_static']
# Add any paths that contain templates here, relative to this directory.
templates_path = ["source/_templates"]

# Output file base name for HTML help builder.
htmlhelp_basename = "project-templatedoc"

# -- Path setup --------------------------------------------------------------
# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

# Insert the current docs directory
sys.path.insert(0, os.path.abspath("."))
# Insert the project root
# sys.path.insert(0, os.path.abspath(".."))
# Insert the package root
sys.path.insert(0, "../kkcalc2/")

os.environ["MPLBACKEND"] = "Agg"  # avoid tkinter import errors on rtfd.io

# -- Project information -----------------------------------------------------

project = "KKCalc2"
copyright = (
    f"2024-{date.today().year}, Matthew Gebert, Benjamin Watts, kkcalc maintainers"
)
author = "Matthew Gebert, Benjamin Watts"

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#
# The short X.Y version.

# version = .__version__
# The full version, including alpha/beta/rc tags.
release = kkcalc2.__version__
version = kkcalc2.__version__
html_title = f"{project} v{version} Manual"

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    # "sphinx.ext.doctest",
    # "sphinx.ext.intersphinx",
    # "sphinx.ext.todo",
    "numpydoc",
    # "sphinx.ext.ifconfig",
    # "sphinx.ext.viewcode",
    # "sphinx.ext.imgmath",
    # "sphinx.ext.napoleon",
    "sphinx_copybutton",
    # "autoapi.extension",
    "matplotlib.sphinxext.plot_directive",
    "sphinx.ext.coverage",
]


# The root toctree document
master_doc = "index"  # NOTE: will be changed to `root_doc` in sphinx 4
root_doc = master_doc

# Setup the API auto-documentation
autosummary_generate = True
autosummary_generate_overwrite = True
# autosummary_imported_members = True
autosummary_ignore_module_all = True

# autoapi_dirs = ["../kkcalc"]
# autoapi_add_toctree_entry = False


numpydoc_xref_param_type = True
numpydoc_xref_ignore = {"optional", "type_without_description", "BadException"}
# Run docstring validation as part of build process
numpydoc_validation_checks = {"all", "GL01", "SA04", "RT03"}
numpydoc_show_class_members = True

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = "sphinx"

# -- Options for LaTeX output ---------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    # 'papersize': 'letterpaper',
    # The font size ('10pt', '11pt' or '12pt').
    # 'pointsize': '10pt',
    # Additional stuff for the LaTeX preamble.
    # 'preamble': '',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (
        "index",
        "kkcalc2.tex",
        "kkcalc2 Documentation",
        "kkcalc2 maintainers",
        "manual",
    ),
]

# -- Intersphinx setup ----------------------------------------------------

# Example configuration for intersphinx: refer to several Python libraries.
# intersphinx_mapping = get_intersphinx_mapping(packages={
#     # "python",
#     "numpy",
# })
