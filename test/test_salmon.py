"""Test importing salmon quantification files."""

from pathlib import Path
from typing import List

import anndata as ad
import pandas as pd
import xarray as xr

from pytximport import tximport
from pytximport.utils import biotype_filters


def test_salmon(
    salmon_file: Path,
    transcript_gene_mapping_human: pd.DataFrame,
) -> None:
    """Test importing a salmon quantification file.

    Args:
        salmon_file (Path): [description]
    """
    for counts_from_abundance in [None, "scaled_tpm", "length_scaled_tpm"]:
        for output_type in ["xarray", "anndata"]:
            result = tximport(
                [salmon_file],
                "salmon",
                transcript_gene_mapping_human,
                output_type=output_type,  # type: ignore
                ignore_transcript_version=True,
                ignore_after_bar=True,
                counts_from_abundance=counts_from_abundance,  # type: ignore
            )

            # check that the result is an xarray dataset
            if output_type == "xarray":
                assert isinstance(result, xr.Dataset)
                counts = result["counts"].data
            else:
                assert isinstance(result, ad.AnnData)
                counts = result.X

            # check that the counts.data are all positive
            assert (counts >= 0).all()


def test_multiple_salmon(
    salmon_multiple_files: List[Path],
    transcript_gene_mapping_mouse: pd.DataFrame,
) -> None:
    """Test importing salmon quantification files.

    Args:
        salmon_multiple_files (Path): [description]
    """
    for counts_from_abundance in [None, "scaled_tpm", "length_scaled_tpm"]:
        for biotype_filter in [None, biotype_filters.GENCODE_PROTEIN_CODING]:
            for existence_optional in [True, False]:
                if existence_optional:
                    salmon_multiple_files_original = salmon_multiple_files.copy()
                    # add a non-existent file
                    salmon_multiple_files = salmon_multiple_files + [
                        salmon_multiple_files.pop(0).with_name("non_existent_file.sf")
                    ]

                result = tximport(
                    salmon_multiple_files,
                    "salmon",
                    transcript_gene_mapping_mouse,
                    ignore_transcript_version=True,
                    ignore_after_bar=True,
                    counts_from_abundance=counts_from_abundance,  # type: ignore
                    biotype_filter=biotype_filter,
                    existence_optional=existence_optional,
                )

                if existence_optional:
                    salmon_multiple_files = salmon_multiple_files_original

            # check that the result is an xarray dataset
            assert isinstance(result, xr.Dataset)

            # check that the counts.data are all positive
            assert (result["counts"].data >= 0).all()
