"""pytximport: Gene-level count estimation from transcript quantification files.

The `pytximport` package provides a Python implementation of the `tximport` R package, which provides an easy-to-use
interface for importing transcript quantification files from various tools (e.g., `salmon`, `kallisto`, `RSEM`) into
Python. The package is designed to work with the output of these tools and provide isoform bias-corrected gene-level
counts for downstream analysis.

The package provides a single function, `tximport`, as the main entry point, as well as a command-line interface
(`pytximport`). Further utility functions (e.g., to create a transcript-to-gene map) are provided through the
`pytximport.utils` module.

`pytximport` can output data as AnnData objects and xarray objects or save the data as a CSV file, enabling seamless
integration with other Python packages such as `PyDESeq2`.

.. note::
    Please consider citing both `pytximport` and the original R implementation of `tximport` when using `pytximport`.
    For `pytximport`, please refer to the README or CITATION.cff file for the appropriate citation information.
    The original `tximport` publication can be found at: https://doi.org/10.12688/f1000research.7563.1

.. warning::
    Though `pytximport` aims to provide the same functionality as `tximport`, there are some differences between the two
    packages. While the same configuration will result in identical output, the configuration options and defaults may
    differ between the two packages. Please refer to the documentation for more information.
"""

from . import definitions, importers, utils
from ._cli import cli
from ._version import __version__
from .core import tximport

# Allow users to import the tximport function as pytximport as well
pytximport = tximport
