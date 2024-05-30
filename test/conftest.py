"""Create example data for pytest."""

import os
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd
import pytest

np.random.seed(1308419541)

FILE_DIR = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "data",
)


@pytest.fixture(scope="session")
def transcript_gene_mapping_human() -> pd.DataFrame:
    """Provide transcript to gene mapping for human samples.

    Returns:
        pd.DataFrame: The transcript to gene mapping.
    """
    return pd.read_table(Path(FILE_DIR) / "gencode.v46.metadata.HGNC.tsv", header=0, sep="\t")


@pytest.fixture(scope="session")
def transcript_gene_mapping_path_mouse() -> Path:
    """Provide the path to the transcript to gene mapping for mouse samples.

    Returns:
        Path: The path to the transcript to gene mapping.
    """
    return Path(FILE_DIR) / "gencode.vM35.metadata.MGI.tsv"


@pytest.fixture(scope="session")
def transcript_gene_mapping_mouse(
    transcript_gene_mapping_path_mouse: Path,
) -> pd.DataFrame:
    """Provide transcript to gene mapping for mouse samples.

    Returns:
        pd.DataFrame: The transcript to gene mapping.
    """
    return pd.read_table(transcript_gene_mapping_path_mouse, header=0, sep="\t")


@pytest.fixture(scope="session")
def kallisto_file() -> Path:
    """Provide the path to a kallisto quantification file.

    The provided path is a directory containing the abundance.h5 file and tests that the importer can handle this.

    Returns:
        Path: The path to the kallisto quantification file.
    """
    return Path(FILE_DIR) / "kallisto"


@pytest.fixture(scope="session")
def kallisto_multiple_files() -> List[Path]:
    """Provide the path to multiple kallisto quantification files.

    These files are .tsv files and therefore also test that part of the kallisto importer.

    Returns:
        List[Path]: The paths to the kallisto quantification files.
    """
    file_paths = []
    for i in range(2):
        file_path = Path(FILE_DIR) / "kallisto" / "multiple" / f"abundance_{i + 1}.tsv"
        file_paths.append(file_path)

    return file_paths


@pytest.fixture(scope="session")
def salmon_file() -> Path:
    """Provide the path to a salmon quantification file."""
    return Path(FILE_DIR) / "salmon" / "quant.sf"


@pytest.fixture(scope="session")
def salmon_multiple_files() -> List[Path]:
    """Create multiple salmon quantification files."""
    file_paths = []
    for i in range(3):
        file_path = Path(FILE_DIR) / "salmon" / "multiple" / f"Sample_{i + 1}.sf"
        file_paths.append(file_path)

    return file_paths
