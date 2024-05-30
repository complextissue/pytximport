"""Configure the Sphinx documentation builder."""

# -- Path setup --------------------------------------------------------------
import os
import sys
from pathlib import Path
from typing import Any, Dict, List

from pygments.styles import get_style_by_name

DOCS_DIR = Path(__file__).parent
sys.path.insert(0, os.path.abspath("../../pytximport"))

# -- Project information -----------------------------------------------------

project = "pytximport"
copyright = "2024, Malte Kuehl"
author = "Malte Kuehl"

# The full version, including alpha/beta/rc tags
release = "0.1.1"

# -- General configuration ---------------------------------------------------

templates_path = ["_templates"]
exclude_patterns = [".DS_Store"]
source_suffix = {
    ".rst": "restructuredtext",
    ".txt": "markdown",
    ".md": "markdown",
}
master_doc = "index"
default_role = "literal"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.napoleon",
    "sphinx.ext.githubpages",
    "sphinx_copybutton",
    "sphinx_design",
    "sphinx_autodoc_typehints",
    "sphinx.ext.inheritance_diagram",
    "nbsphinx",
    "autoapi.extension",
    "myst_parser",
]
myst_enable_extensions = [
    "amsmath",
    "attrs_inline",
    "colon_fence",
    "deflist",
    "dollarmath",
    "fieldlist",
    "html_admonition",
    "html_image",
    "replacements",
    "smartquotes",
    "strikethrough",
    "substitution",
    "tasklist",
]

# don't run the notebooks
nbsphinx_execute = "never"

autoapi_type = "python"
autoapi_add_toctree_entry = False
autoapi_ignore: List[str] = ["_*.py"]
autoapi_options = [
    "members",
    "undoc-members",
    "show-inheritance",
    "show-module-summary",
    "special-members",
    "imported-members",
]
autoapi_dirs = ["../../pytximport"]
autoapi_member_order = "alphabetical"
# autoapi_root = '/api'
autoapi_python_class_content = "init"  # ensures that the __init__ method is also displayed

autodoc_typehints = "description"

napoleon_google_docstring = True
napoleon_numpy_docstring = False
napoleon_include_init_with_doc = True
napoleon_include_special_with_doc = True
napoleon_use_rtype = True
napoleon_use_param = True

intersphinx_mapping = dict(  # noqa: C408
    cycler=("https://matplotlib.org/cycler/", None),
    matplotlib=("https://matplotlib.org/", None),
    numpy=("https://docs.scipy.org/doc/numpy/", None),
    pandas=("https://pandas.pydata.org/pandas-docs/stable/", None),
    pytest=("https://docs.pytest.org/en/latest/", None),
    python=("https://docs.python.org/3", None),
    scipy=("https://docs.scipy.org/doc/scipy/reference/", None),
    seaborn=("https://seaborn.pydata.org/", None),
    sklearn=("https://scikit-learn.org/stable/", None),
)

# -- Options for HTML output -------------------------------------------------

html_theme = "furo"

html_title = "pytximport"

html_static_path = ["_static"]
html_logo = "_static/pytximport_logo.png"

html_css_files = [
    "css/custom.css",
]

html_context: Dict[str, Any] = {
    "display_github": True,
    "github_user": "complextissue",
    "github_repo": "pytximport",
    "github_version": "main",
    "conf_py_path": "/docs/source/",
    "github_button": True,
    "show_powered_by": False,
}
html_show_sphinx = False

pygments_style = "monokai"

plot_include_source = True
plot_formats = [("svg", 90)]
plot_html_show_formats = False
plot_html_show_source_link = False
plot_working_directory = DOCS_DIR.parent
