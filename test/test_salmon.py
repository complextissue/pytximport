"""Test importing salmon quantification files."""

from pathlib import Path
from typing import List

import anndata as ad
import pandas as pd
import xarray as xr
from summarizedexperiment import SummarizedExperiment

from pytximport import tximport
from pytximport.utils import biotype_filters


def test_salmon(
    salmon_file: Path,
    transcript_gene_mapping_human: pd.DataFrame,
) -> None:
    """Test importing a salmon quantification file.

    Args:
        salmon_file (Path): Path to the salmon quantification file.
        transcript_gene_mapping_human (pd.DataFrame): Transcript to gene mapping.
    """
    for counts_from_abundance in [None, "scaled_tpm", "length_scaled_tpm"]:
        for output_type in ["xarray", "anndata", "summarizedexperiment"]:
            result = tximport(
                [salmon_file],
                "salmon",
                transcript_gene_mapping_human,
                output_type=output_type,  # type: ignore
                ignore_transcript_version=True,
                ignore_after_bar=True,
                counts_from_abundance=counts_from_abundance,  # type: ignore
            )

            # Check that the result is an xarray dataset
            if output_type == "xarray":
                assert isinstance(result, xr.Dataset)
                counts = result["counts"].data
            elif output_type == "summarizedexperiment":
                assert isinstance(result, SummarizedExperiment)
                counts = result.assays["counts"]
            else:
                assert isinstance(result, ad.AnnData)
                counts = result.X

            # Check that the counts.data are all positive
            assert (counts >= 0).all()


def test_salmon_gzip(
    salmon_file_gzip: Path,
    transcript_gene_mapping_human: pd.DataFrame,
) -> None:
    """Test importing a gzipped salmon quantification file.

    Args:
        salmon_file_gzip (Path): Path to the gzipped salmon quantification file.
        transcript_gene_mapping_human (pd.DataFrame): Transcript to gene mapping.
    """
    result = tximport(
        [salmon_file_gzip],
        "salmon",
        transcript_gene_mapping_human,
        ignore_transcript_version=True,
        ignore_after_bar=True,
    )

    # Check that the result is an AnnData object
    assert isinstance(result, ad.AnnData)

    # Check that the counts.data are all positive
    assert (result.X >= 0).all()


def test_multiple_salmon(
    salmon_multiple_files: List[Path],
    transcript_gene_mapping_mouse: pd.DataFrame,
) -> None:
    """Test importing salmon quantification files.

    Args:
        salmon_multiple_files (Path): List of paths to the salmon quantification files.
        transcript_gene_mapping_mouse (pd.DataFrame): Transcript to gene mapping.
    """
    for counts_from_abundance in [None, "scaled_tpm", "length_scaled_tpm"]:
        for biotype_filter in [None, biotype_filters.GENCODE_PROTEIN_CODING]:
            for existence_optional in [True, False]:
                if existence_optional:
                    salmon_multiple_files_original = salmon_multiple_files.copy()
                    # Add a non-existent file
                    salmon_multiple_files = salmon_multiple_files + [
                        salmon_multiple_files.pop(0).with_name("non_existent_file.sf")
                    ]

                result = tximport(
                    salmon_multiple_files,
                    "salmon",
                    transcript_gene_mapping_mouse,
                    ignore_transcript_version=True,
                    ignore_after_bar=True,
                    output_type="xarray",
                    counts_from_abundance=counts_from_abundance,  # type: ignore
                    biotype_filter=biotype_filter,
                    existence_optional=existence_optional,
                )

                if existence_optional:
                    salmon_multiple_files = salmon_multiple_files_original

            # Check that the result is an xarray dataset
            assert isinstance(result, xr.Dataset)

            # Check that the counts.data are all positive
            assert (result["counts"].data >= 0).all()
