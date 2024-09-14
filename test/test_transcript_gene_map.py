"""Test importing salmon quantification files."""

from pathlib import Path

import pandas as pd

from pytximport.utils import (
    create_transcript_gene_map,
    create_transcript_gene_map_from_annotation,
)


def test_transcript_to_gene_map() -> None:
    """Test creating a transcript to gene map."""
    for species in ["human", "mouse"]:
        df_transcript_to_gene = create_transcript_gene_map(
            species=species,
            host="http://www.ensembl.org",
            target_field="external_gene_name",
        )

        assert isinstance(df_transcript_to_gene, pd.DataFrame), "The output is not a DataFrame."
        assert df_transcript_to_gene.shape[1] == 2, "The output has the wrong number of columns."
        assert "transcript_id" in df_transcript_to_gene.columns, "The output is missing the `transcript_id` column."
        assert "gene_id" in df_transcript_to_gene.columns, "The output is missing the `gene_id` column."


def test_transcript_to_gene_map_from_gtf_annotation(
    gtf_annotation_file: Path,
) -> None:
    """Test creating a transcript to gene map from a GTF annotation file."""
    for use_gene_name in [True, False]:
        for keep_biotype in [True, False]:
            df_transcript_to_gene = create_transcript_gene_map_from_annotation(
                gtf_annotation_file,
                target_field=("gene_id" if not use_gene_name else "gene_name"),
                keep_biotype=keep_biotype,
            )

            assert isinstance(df_transcript_to_gene, pd.DataFrame), "The output is not a DataFrame."

            if use_gene_name:
                # Check that not all gene ids start with ENSG
                assert (
                    not df_transcript_to_gene["gene_id"].str.startswith("ENSG").all()
                ), "All gene ids start with ENSG."

            if keep_biotype:
                assert df_transcript_to_gene.shape[1] == 3, "The output has the wrong number of columns."
            else:
                assert df_transcript_to_gene.shape[1] == 2, "The output has the wrong number of columns."

            df_transcript_to_gene_reference = pd.read_csv(
                gtf_annotation_file.parent / f"transcript_to_gene_map{'_gene_name' if use_gene_name else ''}.csv",
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

    # Test using the transcript name as the source field
    df_transcript_to_gene = create_transcript_gene_map_from_annotation(
        gtf_annotation_file,
        source_field="transcript_name",
        target_field="gene_id",
    )

    assert isinstance(df_transcript_to_gene, pd.DataFrame), "The output is not a DataFrame."
    assert df_transcript_to_gene.shape[1] == 2, "The output has the wrong number of columns."
