"""Test importing salmon quantification files."""

from pathlib import Path
from typing import List

import anndata as ad
import pandas as pd

from pytximport.utils import (
    create_transcript_to_gene_map,
    create_transcript_to_gene_map_from_gtf_annotation,
)


def test_transcript_to_gene_map() -> None:
    """Test creating a transcript to gene map."""
    df_transcript_to_gene = create_transcript_to_gene_map(
        species="human",
        host="http://www.ensembl.org",
        field="external_gene_name",
    )

    assert isinstance(df_transcript_to_gene, pd.DataFrame), "The output is not a DataFrame."
    assert df_transcript_to_gene.shape[1] == 2, "The output has the wrong number of columns."
    assert "transcript_id" in df_transcript_to_gene.columns, "The output is missing the `transcript_id` column."
    assert "gene_id" in df_transcript_to_gene.columns, "The output is missing the `gene_id` column."


def test_transcript_to_gene_map_from_gtf_annotation(
    gtf_annotation_file: Path,
) -> None:
    """Test creating a transcript to gene map from a GTF annotation file."""
    for keep_gene_name in [True, False]:
        for keep_biotype in [True, False]:
            df_transcript_to_gene = create_transcript_to_gene_map_from_gtf_annotation(
                gtf_annotation_file,
                field="gene_id",
                keep_gene_name=keep_gene_name,
                keep_biotype=keep_biotype,
            )

            assert isinstance(df_transcript_to_gene, pd.DataFrame), "The output is not a DataFrame."

            if keep_gene_name and keep_biotype:
                assert df_transcript_to_gene.shape[1] == 4, "The output has the wrong number of columns."
            elif keep_gene_name or keep_biotype:
                assert df_transcript_to_gene.shape[1] == 3, "The output has the wrong number of columns."
            else:
                assert df_transcript_to_gene.shape[1] == 2, "The output has the wrong number of columns."

            df_transcript_to_gene_reference = pd.read_csv(
                gtf_annotation_file.parent / f"transcript_to_gene_map{'_gene_name' if keep_gene_name else ''}.csv",
                header=0,
            ).reset_index(drop=True)

            pd.testing.assert_frame_equal(
                (
                    df_transcript_to_gene.drop(columns=["gene_biotype"])
                    if "gene_biotype" in df_transcript_to_gene.columns
                    else df_transcript_to_gene
                ),
                df_transcript_to_gene_reference,
            )
