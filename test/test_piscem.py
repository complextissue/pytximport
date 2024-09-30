"""Test importing salmon quantification files."""

from pathlib import Path
from typing import List

import anndata as ad
import pandas as pd
import xarray as xr

from pytximport import tximport
from pytximport.utils import biotype_filters


def test_piscem(
    piscem_file: Path,
    transcript_gene_mapping_human: pd.DataFrame,
) -> None:
    """Test importing a piscem quantification file.

    Args:
        piscem_file (Path): Path to the piscem quantification file.
        transcript_gene_mapping_human (pd.DataFrame): Transcript to gene mapping.
    """
    result = tximport(
        [piscem_file],
        "piscem",
        transcript_gene_mapping_human,
        output_type="anndata",
        ignore_transcript_version=True,
        ignore_after_bar=True,
        counts_from_abundance=None,  # type: ignore
    )

    # Check that the result is an AnnData object
    assert isinstance(result, ad.AnnData)

    # Check that the counts.data are all positive
    assert (result.X >= 0).all()
