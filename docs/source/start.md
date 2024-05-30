# Start

[![Version](https://img.shields.io/pypi/v/pytximport)](https://pypi.org/project/pytximport/)
[![Code Style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)](https://github.com/pre-commit/pre-commit)
[![Python Version Required](https://img.shields.io/pypi/pyversions/pytximport)](https://pypi.org/project/pytximport/)
[![License](https://img.shields.io/pypi/l/pytximport)](https://github.com/complextissue/pytximport)

`pytximport` is a Python package for fast gene count estimation based on transcript abundance, inspired by the `tximport` R package.

## Development status

`pytximport` is alpha version software. While it should work for most use cases we regularly compare outputs against the R implementation, expect breaking changes and bugs. If you encounter any problems, please open a GitHub issue. If you are a Python developer, we welcome pull requests implementing missing features, adding more extensive unit tests and bug fixes.

## Installation

```bash
pip install pytximport
```

## Quick Start

You can either use it from the command line:

```bash
pytximport -i ./sample_1.sf -i ./sample_2.sf -t salmon -m ./tx2gene_map.tsv -o ./output_counts.csv
```

Common options are:
- `-i`: The input files.
- `-t`: The input type, e.g., `salmon`, `kallisto` or `tsv`.
- `-m`: The map to match transcript ids to their gene ids. Expected column names are `transcript_id` and `gene_id`.
- `-o`: The output path.
- `-c`: The count transform to apply. Leave out for none, other options include `scaled_tpm` and `length_scaled_tpm`.
- `-id`: The column name containing the transcript ids, in case it differs from the typical naming standards for the configured input file type.
- `-counts`: The column name containing the transcript counts, in case it differs from the typical naming standards for the configured input file type.
- `-length`: The column name containing the transcript lenghts, in case it differs from the typical naming standards for the configured input file type.
- `-tpm`: The column name containing the transcript abundance, in case it differs from the typical naming standards for the configured input file type.

Or import the `tximport` function in your Python files:

```python
from pytximport import tximport
results = tximport(file_paths, "salmon", transcript_gene_mapping)
```

## Motivation

The `tximport` package has become a main stay in the bulk RNA sequencing community and has been used in hundreds of scientific publications. However, its accessibility has remained limited since it requires the R programming language and cannot be used from within Python scripts or the command line. Other tools of the bulk RNA sequencing analysis stack, like `DESeq2` (in the form of `PyDESeq2`), `decoupler`, `liana` and others all have Python versions. Additionally, pseudoalignment tools like `salmon` and `kallisto` can be installed via `conda` and can be used from the command line.
`tximport` thus constitutes the missing link in many common analysis workflows. `pytximport` fills this gap and allows these workflows to be entirely done in Python, which is preinstalled on most development machines, and from the command line.

## Citation

Please cite both the original publication as well as this Python implementation:
- Charlotte Soneson, Michael I. Love, Mark D. Robinson. Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences, F1000Research, 4:1521, December 2015. doi: 10.12688/f1000research.7563.1
- Kuehl, M., & Puelles, V. (2024). pytximport: Fast gene count estimation from transcript quantification files in Python (Version 0.1.1) [Computer software]. https://github.com/complextissue/pytximport

## Differences

Generally, outputs from `pytximport` correspond to the outputs from `tximport` within the accuracy allowed by multiple floating point operations and small implementation differences in its dependencies when using the same configuration. If you observe larger discrepancies, please open an issue.

While the outputs are roughly identical for the same configuration, there remain some differences between the packages:
- `pytximport` can be used from the command line.
- `pytximport` supports `AnnData` format outputs (set `output_type` to `anndata`), enabling seamless integration with the `scverse`.
- `pytximport` currently does not support inferential replicates. If these are valuable to your workflow, we appreciate pull requests to add support.
- `pytximport` currently does not support gene-level inputs. If these are valuable to your workflow, we appreciate pull requests to add support.
- Argument order and argument defaults may differ between the implementations.
- Additional features:
    - When `ignore_transcript_version` is set, the transcript version will not only be scrapped from the quantization file but also from the provided transcript to gene mapping.
    - When `biotype_filter` is set, all transcripts that do not contain any of the provided biotypes will be removed prior to all other steps.
    - When `save_path` is configured, a count matrix will be saved as a .csv file.

## License

The software is provided under the GNU Public License version 3. Please consult `LICENSE.md` for further information.
