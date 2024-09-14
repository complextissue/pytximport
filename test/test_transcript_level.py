"""Test importing salmon quantification files."""

from pathlib import Path
from typing import List

import anndata as ad
import pandas as pd

from pytximport import tximport
from pytximport.utils import biotype_filters, replace_transcript_ids_with_names


def test_salmon_transcript_level(
    salmon_file: Path,
    transcript_name_mapping_human: pd.DataFrame,
    transcript_name_mapping_human_path: Path,
) -> None:
    """Test importing a salmon quantification file at the transcript level.

    Args:
        salmon_file (Path): Path to the salmon quantification file.
        transcript_name_mapping_human (pd.DataFrame): The mapping of transcript ids to transcript names.
        transcript_name_mapping_human_path (Path): The path to the transcript name mapping.
    """
    result = tximport(
        [salmon_file],
        "salmon",
        return_transcript_data=True,
        output_type="anndata",
        ignore_transcript_version=True,
        ignore_after_bar=True,
        counts_from_abundance="scaled_tpm",
    )

    assert isinstance(result, ad.AnnData)
    counts = result.X

    # Check that the counts.data are all positive
    assert (counts >= 0).all()

    # replace the transcript ids with the transcript names
    result_df = replace_transcript_ids_with_names(result, transcript_name_mapping_human)

    # Check that the result is an AnnData object
    assert isinstance(result_df, ad.AnnData)

    # Check that the var names don't start with ENST
    assert result_df.var_names.str.startswith("ENST").sum() == 0

    # Check that the var names are not nan or the string "nan" or empty
    assert (
        result_df.var_names.isna().sum() == 0
        and (result_df.var_names == "nan").sum() == 0
        and (result_df.var_names == "").sum() == 0
    )

    # Check that it also works with a path to the transcript name mapping
    result_path = replace_transcript_ids_with_names(result, transcript_name_mapping_human_path)

    # Check that the results are the same
    pd.testing.assert_frame_equal(result_df.var_names.to_frame(), result_path.var_names.to_frame())


def test_dtu_scaled_tpm(
    salmon_file: Path,
    salmon_multiple_files: List[Path],
    transcript_gene_mapping_human: pd.DataFrame,
    transcript_gene_mapping_mouse: pd.DataFrame,
) -> None:
    """Test importing a salmon quantification file with DTU scaled TPM.

    Args:
        salmon_file (Path): Path to the salmon quantification file.
        transcript_gene_mapping_human (pd.DataFrame): The mapping of transcript ids to gene ids.
    """
    result = tximport(
        [salmon_file],
        "salmon",
        transcript_gene_mapping_human,
        return_transcript_data=True,
        output_type="anndata",
        ignore_transcript_version=True,
        ignore_after_bar=True,
        counts_from_abundance="dtu_scaled_tpm",
    )

    assert isinstance(result, ad.AnnData)
    assert (result.X >= 0).all()

    result = tximport(
        salmon_multiple_files,
        "salmon",
        transcript_gene_mapping_mouse,
        return_transcript_data=True,
        output_type="anndata",
        ignore_transcript_version=True,
        ignore_after_bar=True,
        counts_from_abundance="dtu_scaled_tpm",
        # dtu_scaled_tpm requires all transcript ids to have a gene_id, so we filter out the non-protein-coding
        biotype_filter=biotype_filters.GENCODE_PROTEIN_CODING,
    )

    assert isinstance(result, ad.AnnData)
    assert (result.X >= 0).all()
