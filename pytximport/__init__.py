"""pytximport: Convert transcript-level expression to gene-level expression.

The `pytximport` package provides a Python implementation of the `tximport` R package, which converts transcript-level
expression to gene-level count estimates.

.. note::
    Please also cite the original `tximport` R package when using `pytximport`.
    DOI: https://doi.org/10.12688/f1000research.7563.1

.. warning::
    The `pytximport` package is still in development and currently does not make use of the inferential replicates,
    unlike the `tximport` R package. Default parameters, argument names and return types may also differ between the
    implementations.
"""

from . import definitions, importers, utils
from ._cli import cli
from ._version import __version__
from .core import tximport

# allow users to import the tximport function as pytximport as well
pytximport = tximport
