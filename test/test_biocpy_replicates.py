"""Test importing salmon quantification files."""

from pathlib import Path
from typing import List

from summarizedexperiment import SummarizedExperiment

from pytximport import tximport


def test_biocpy_inferential_replicates(
    fabry_disease_files: List[Path],
) -> None:
    """Test importing quantification files with inferential replicates with SummarizedExperiment output.

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
        output_type="summarizedexperiment",
        inferential_replicates=True,
        inferential_replicate_variance=True,
        inferential_replicate_transformer=None,  # Do not recalculate counts but include the inferential replicates
    )

    assert isinstance(result, SummarizedExperiment), "The result is not an SummarizedExperiment object."

    # Check that variance is in the metadata
    assert "variance" in result.get_metadata()

    # Check that inferential replicates is in the metadata
    assert "inferential_replicates" in result.get_metadata()
