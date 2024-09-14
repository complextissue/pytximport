"""Test importing salmon quantification files."""

from pathlib import Path
from typing import List

import anndata as ad

from pytximport import tximport


def test_anndata_replicates(
    fabry_disease_files: List[Path],
) -> None:
    """Test importing quantification files with inferential replicates with AnnData output.

    Args:
        fabry_disease_files (List[Path]): The paths to the quantification files.
    """
    fabry_directory = fabry_disease_files[0].parent

    result = tximport(
        fabry_disease_files,
        "salmon",
        fabry_directory / "transcript_gene_mapping_human.tsv",
        ignore_transcript_version=True,
        ignore_after_bar=True,
        output_type="anndata",
        inferential_replicates=True,
        inferential_replicate_variance=True,
        inferential_replicate_transformer=None,  # Do not recalculate counts but include the inferential replicates
    )

    assert isinstance(result, ad.AnnData), "The result is not an AnnData object."

    # Check that variance is in the obsm
    assert "variance" in result.obsm.keys()

    # Check that inferential replicates is in the uns
    assert "inferential_replicates" in result.uns.keys()
