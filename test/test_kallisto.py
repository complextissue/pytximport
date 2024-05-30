"""Test importing kallisto quantification files."""

from pathlib import Path
from typing import List

import pandas as pd
import xarray as xr

from pytximport import tximport


def test_kallisto(
    kallisto_file: Path,
    transcript_gene_mapping_human: pd.DataFrame,
) -> None:
    """Test importing a kallisto quantification file.

    Args:
        kallisto_file (Path): [description]
    """
    for counts_from_abundance in [None, "scaled_tpm", "length_scaled_tpm"]:
        result = tximport(
            [kallisto_file],
            "kallisto",
            transcript_gene_mapping_human,
            ignore_transcript_version=True,
            ignore_after_bar=True,
            counts_from_abundance=counts_from_abundance,  # type: ignore
        )

        # check that the result is an xarray dataset
        assert isinstance(result, xr.Dataset)

        # check that the counts.data are all positive
        assert (result["counts"].data >= 0).all()


def test_multiple_kallisto(
    kallisto_multiple_files: List[Path],
    transcript_gene_mapping_human: pd.DataFrame,
) -> None:
    """Test importing kallisto quantification files.

    Args:
        kallisto_multiple_files (Path): [description]
    """
    for counts_from_abundance in [None, "scaled_tpm", "length_scaled_tpm"]:
        for existence_optional in [True, False]:
            if existence_optional:
                kallisto_multiple_files_original = kallisto_multiple_files.copy()
                # add a non-existent file
                kallisto_multiple_files = kallisto_multiple_files + [
                    kallisto_multiple_files.pop(0).with_name("non_existent_file.tsv")
                ]

            result = tximport(
                kallisto_multiple_files,
                "kallisto",
                transcript_gene_mapping_human,
                id_column="target_id",
                counts_column="est_counts",
                length_column="eff_length",
                abundance_column="tpm",
                ignore_transcript_version=True,
                ignore_after_bar=True,
                counts_from_abundance=counts_from_abundance,  # type: ignore
                existence_optional=existence_optional,
            )

            if existence_optional:
                kallisto_multiple_files = kallisto_multiple_files_original

        # check that the result is an xarray dataset
        assert isinstance(result, xr.Dataset)

        # check that the counts.data are all positive
        assert (result["counts"].data >= 0).all()
