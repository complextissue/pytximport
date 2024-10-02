# pytximport

<hr />

[![Version](https://img.shields.io/pypi/v/pytximport)](https://pypi.org/project/pytximport/)
[![License](https://img.shields.io/pypi/l/pytximport)](https://github.com/complextissue/pytximport)
![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/complextissue/pytximport/ci.yml)
[![Documentation Status](https://readthedocs.org/projects/pytximport/badge/?version=latest)](https://pytximport.readthedocs.io/en/latest/?badge=latest)
[![Codecov](https://codecov.io/gh/complextissue/pytximport/graph/badge.svg?token=M9JEHJVXYI)](https://codecov.io/gh/complextissue/pytximport)
[![Install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/pytximport/README.html)
![PyPI - Downloads](https://img.shields.io/pypi/dm/pytximport)
[![Python Version Required](https://img.shields.io/pypi/pyversions/pytximport)](https://pypi.org/project/pytximport/)
[![Code Style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)](https://github.com/pre-commit/pre-commit)

`pytximport` is a Python package for efficient (gene-)count estimation from transcript quantification files produced by pseudoalignment/quasi-mapping tools such as `salmon`, `kallisto`, `rsem` and others. `pytximport` is a port of the popular [tximport Bioconductor R package](https://bioconductor.org/packages/release/bioc/html/tximport.html).

## Installation

The recommended way to install `pytximport` is through Bioconda:

```bash
mamba install -c bioconda pytximport
```

`pytximport` can also be installed via pip:

```bash
pip install pytximport
```

While not required, we recommend users also install `pyarrow` for faster import of tab-separated value-based quantification files:

```bash
mamba install -c conda-forge pyarrow-core
```

or:

```bash
pip install pyarrow
```

## Quick Start

You can either import the `tximport` function in your Python files:

```python
from pytximport import tximport
from pytximport.utils import create_transcript_gene_map

transcript_gene_map = create_transcript_gene_map(species="human")

results = tximport(
    file_paths,
    data_type="salmon",
    transcript_gene_map=transcript_gene_map,
)
```

Or use it from the command line:

```bash
pytximport -i ./sample_1.sf -i ./sample_2.sf -t salmon -m ./tx2gene_map.tsv -o ./output_counts.csv
```

Common options are:

- `-i`: The path to an quantification file. To provide multiple input files, use `-i input1.sf -i input2.sf ...`.
- `-t`: The type of quantification file, e.g. `salmon`, `kallisto` and others.
- `-m`: The path to the transcript to gene map. Either a tab-separated (.tsv) or comma-separated (.csv) file. Expected column names are `transcript_id` and `gene_id`.
- `-o`: The output path to save the resulting counts to.
- `-of`: The format of the output file. Either `csv` or `h5ad`.
- `-ow`: Provide this flag to overwrite an existing file at the output path.
- `-c`: The method to calculate the counts from the abundance. Leave empty to use counts. For differential gene expression analysis, we recommend using `length_scaled_tpm`. For differential transcript expression analysis, we recommend using `scaled_tpm`. For differential isoform usage analysis, we recommend using `dtu_scaled_tpm`.
- `-ir`: Provide this flag to make use of inferential replicates. Will use the median of the inferential replicates.
- `-gl`: Provide this flag when importing gene-level counts from RSEM files.
- `-tx`: Provide this flag to return transcript-level instead of gene-summarized data. Incompatible with gene-level input and `counts_from_abundance=length_scaled_tpm`.
- `--help`: Display all configuration options.

## Documentation

Detailled documentation is made available at: [https://pytximport.readthedocs.io](https://pytximport.readthedocs.io/en/latest/start.html).

## Development status

`pytximport` is still in development and has not yet reached version 1.0.0 in the [SemVer](https://semver.org/) versioning scheme. While it should work for almost all use cases and we regularly compare outputs against the R implementation, breaking changes between minor versions may occur. If you encounter any problems, please open a GitHub issue. If you are a Python developer, we welcome pull requests implementing missing features, adding more extensive unit tests and bug fixes.

## Motivation

The `tximport` package has become a main stay in the bulk RNA sequencing community and has been used in hundreds of scientific publications. However, its accessibility has remained limited since it requires the R programming language and cannot be used from within Python scripts or the command line. Other tools of the bulk RNA sequencing analysis stack, like `DESeq2` (in the form of `PyDESeq2`), `decoupler`, `liana` and others all have Python versions. Additionally, pseudoalignment tools like `salmon` and `kallisto` can be installed via `conda` and can be used from the command line.
`tximport` thus constitutes the missing link in many common analysis workflows. `pytximport` fills this gap and allows these workflows to be entirely done in Python, which is preinstalled on most development machines, and from the command line.

## Citation

Please cite both the original publication as well as this Python implementation:

- Kuehl, M., & Puelles, V. (2024). pytximport: Gene count estimation from transcript quantification files in Python (Version 0.10.0) [Computer software]. https://github.com/complextissue/pytximport
- Charlotte Soneson, Michael I. Love, Mark D. Robinson. Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences, F1000Research, 4:1521, December 2015. doi: 10.12688/f1000research.7563.1

## License

The software is provided under the GNU General Public License version 3. Please consult `LICENSE` for further information.

## Differences

Generally, outputs from `pytximport` correspond to the outputs from `tximport` within the accuracy allowed by multiple floating point operations and small implementation differences in its dependencies when using the same configuration. If you observe larger discrepancies, please open an issue.

While the outputs are identical within floating point tolerance for the same configuration, there remain some differences between the packages:

Features unique to `pytximport`:
- Generating transcript-to-gene maps, either from a BioMart server or an `annotation.gtf` file. Use `create_transcript_gene_map` or `create_transcript_gene_map_from_annotation` from `pytximport.utils`.
- Command line interface. Type `pytximport --help` into your terminal to explore all options.
- `AnnData`-support, enabling seamless integration with the `scverse`.
- `SummarizedExperiment`-support to represent outputs in familiar Bioconductor data structures available through the [BiocPy](https://github.com/biocpy) ecosystem.
- Saving outputs directly to file (use the `output_path` argument).
- Removing transcript versions from **both** the quantification files and the transcript-to-gene map when `ignore_transcript_version` is provided.
- Post-hoc biotype-filtering. Set `biotype_filter` to a whitelist of possible biotypes contained within the bar-separated values of your transcript ids.

Features unique to `tximport`:
- Alevin single-cell RNA-seq data support

Argument order and argument defaults may differ between the implementations.

## Contributing

Contributions are welcome. Contributors are asked to follow the Contributor Covenant Code of Conduct.

To set up `pytximport` for development on your machine, we recommend to git clone the dev branch:

```bash
git clone --depth 1 -b dev https://github.com/complextissue/pytximport.git
cd pytximport
pyenv local 3.9
make create-venv
source .venv/source/activate
make install-dev
```

Since `pytximport` is linted and formatted, the repository contains a list of recommended VS Code extensions in `.vscode/extensions.json`. If you are using a different editor, please make sure to set up your environment to use the same linters and formatters.

For new features and non-obvious bug fixes, we kindly ask that you create a GitHub issue before submitting a PR.

## Running the tests locally

Please follow the steps described in the "Contributing" section. Once you have setup your development environment, you can run the unit tests locally:

```bash
make coverage-report
```

## Building the documentation locally

The documentation can be build locally by navigating to the `docs` folder and running: `make html`.
This requires that the development requirements of the package as well as the package itself have been installed in the same virtual environment and that `pandoc` has been added, e.g. by running `brew install pandoc` on macOS operating systems.

## Data sources

The quantification files used for the unit tests are partly adopted from [tximportData](https://doi.org/doi:10.18129/B9.bioc.tximportData) which in turn used a subsample of the GEUVADIS data:
Lappalainen, T., Sammeth, M., Friedländer, M. R., ‘t Hoen, P. A., Monlong, J., Rivas, M. A., ... & Dermitzakis, E. T. (2013). Transcriptome and genome sequencing uncovers functional variation in humans. Nature, 501(7468), 506-511.

Other test and example files, such as those used in the vignette, are based on the following work:
Braun, F., Abed, A., Sellung, D., Rogg, M., Woidy, M., Eikrem, O., ... & Huber, T. B. (2023). Accumulation of α-synuclein mediates podocyte injury in Fabry nephropathy. The Journal of clinical investigation, 133(11).
