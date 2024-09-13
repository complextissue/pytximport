"""Test importing RSEM quantification files."""

from pathlib import Path
from typing import List

import pandas as pd
import pytest
import xarray as xr

from pytximport import tximport


def test_rsem(
    rsem_files: List[Path],
    transcript_gene_mapping_human: pd.DataFrame,
) -> None:
    """Test importing RSEM quantification files.

    Args:
        rsem_files (Path): List of paths to the RSEM quantification files.
        transcript_gene_mapping_human (pd.DataFrame): The transcript to gene mapping.
    """
    for gene_level in [True, False]:
        for counts_from_abundance in [None, "length_scaled_tpm"]:
            for existence_optional in [True, False]:
                if existence_optional:
                    rsem_files_original = rsem_files.copy()
                    # Add a non-existent file
                    rsem_files = rsem_files + [rsem_files[-1].with_name("non_existent_file.genes.results.gz")]

                if gene_level and counts_from_abundance is not None:
                    with pytest.raises(ValueError):
                        result = tximport(
                            rsem_files,
                            "rsem",
                            transcript_gene_mapping_human,
                            counts_from_abundance=counts_from_abundance,
                            gene_level=gene_level,
                            ignore_transcript_version=True,
                            ignore_after_bar=True,
                            output_type="xarray",
                            existence_optional=existence_optional,
                        )
                else:
                    result = tximport(
                        rsem_files,
                        "rsem",
                        transcript_gene_mapping_human,
                        counts_from_abundance=counts_from_abundance,
                        gene_level=gene_level,
                        ignore_transcript_version=True,
                        ignore_after_bar=True,
                        output_type="xarray",
                        existence_optional=existence_optional,
                    )

                    # Check that the result is an xarray dataset
                    assert isinstance(result, xr.Dataset)

                    # Check that the counts.data are all positive
                    assert (result["counts"].data >= 0).all()

                    # Check that the gene_id coordinates do not contain any "."
                    assert not any("." in gene_id for gene_id in result.coords["gene_id"].values)

                if existence_optional:
                    rsem_files = rsem_files_original

    # Test gene level without transcript_gene_map
    result = tximport(
        rsem_files,
        "rsem",
        transcript_gene_map=None,
        counts_from_abundance=None,
        gene_level=True,
        output_type="xarray",
    )

    # Check that the result is an xarray dataset
    assert isinstance(result, xr.Dataset)

    # Check that the counts.data are all positive
    assert (result["counts"].data >= 0).all()
