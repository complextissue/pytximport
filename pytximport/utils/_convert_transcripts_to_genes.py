from logging import log, warning
from typing import List, Literal, Optional, Union

import numpy as np
import pandas as pd
import xarray as xr

from ._convert_abundance_to_counts import convert_abundance_to_counts
from ._replace_missing_average_transcript_length import (
    replace_missing_average_transcript_length,
)


def convert_transcripts_to_genes(
    transcript_data: xr.Dataset,
    transcript_gene_map: pd.DataFrame,
    ignore_after_bar: bool = True,
    ignore_transcript_version: bool = True,
    biotype_filter: Optional[List[str]] = None,
    counts_from_abundance: Optional[Literal["scaled_tpm", "length_scaled_tpm"]] = None,
) -> xr.Dataset:
    """Convert transcript-level expression to gene-level expression.

    Args:
        transcript_data (xr.Dataset): The transcript-level expression data from multiple samples.
        transcript_gene_map (pd.DataFrame): The mapping from transcripts to genes. Contains two columns: `transcript_id`
            and `gene_id`.
        ignore_after_bar (bool, optional): Whether to split the transcript id after the bar character (`|`).
            Defaults to True.
        ignore_transcript_version (bool, optional): Whether to ignore the transcript version. Defaults to True.
        counts_from_abundance (Optional[Literal["scaled_tpm", "length_scaled_tpm"]], optional): The type of counts to
            convert to. Defaults to "length_scaled_tpm".

    Returns:
        xr.Dataset: The gene-level expression data from multiple samples.
    """
    transcript_ids: Union[np.ndarray, List[str]] = transcript_data.coords["transcript_id"].values
    transcript_ids = transcript_data.coords["transcript_id"].values

    if biotype_filter is not None:
        # this only works if the transcript_id contains the biotype as one of the bar-separated fields
        assert any(
            "|" in transcript_id for transcript_id in transcript_ids
        ), "The transcript_id does not contain the biotype."

        transcript_id_fields = [transcript_id.split("|") for transcript_id in transcript_ids]

        # iterate through the biotypes and mark a transcript to be kept if a biotype is a match
        transcript_keep_boolean = np.zeros(len(transcript_id_fields))

        for biotype in biotype_filter:
            transcript_keep_boolean = np.logical_or(
                transcript_keep_boolean,
                [(biotype in transcript_id_field) for transcript_id_field in transcript_id_fields],
            )

        # check that at least one transcript is protein-coding
        assert any(transcript_keep_boolean), "No transcripts with the desired biotype are present in the data."

        transcript_data = transcript_data.isel(
            transcript_id=transcript_keep_boolean,
            drop=True,
        )
        log(
            25,
            f"Removed {len(transcript_keep_boolean) - sum(transcript_keep_boolean)} transcripts with other biotypes.",
        )
        transcript_ids = transcript_data.coords["transcript_id"].values

    if ignore_after_bar:
        # ignore the part of the transcript ID after the bar
        transcript_ids = [transcript_id.split("|")[0] for transcript_id in transcript_ids]
        transcript_data.coords["transcript_id"] = transcript_ids

    if ignore_transcript_version:
        # ignore the transcript version in both the data and the transcript gene map
        transcript_gene_map["transcript_id"] = transcript_gene_map["transcript_id"].str.split(".").str[0]
        transcript_ids = [transcript_id.split(".")[0] for transcript_id in transcript_ids]
        transcript_data.coords["transcript_id"] = transcript_ids

    unique_transcripts = list(set(transcript_ids))

    # avoid duplicates in the mapping
    transcripts_duplicated = transcript_gene_map["transcript_id"].duplicated()
    assert not any(transcripts_duplicated), "The mapping contains duplicates."

    # check that at least one transcript is present in the mapping
    assert any(transcript_gene_map["transcript_id"].isin(unique_transcripts)), (
        "No transcripts are present in the mapping. "
        "Please make sure you are using the correct mapping and that the transcript IDs match the mapping. "
        "Adjust the `ignore_after_bar` and `ignore_transcript_version` parameters if necessary."
    )

    # check whether there are any missing transcripts, and if so, warn the user and remove them
    if not set(unique_transcripts).issubset(set(transcript_gene_map["transcript_id"])):
        warning(
            "Not all transcripts are present in the mapping."
            + f" {len(set(unique_transcripts) - set(transcript_gene_map['transcript_id']))}"
            + f" out of {len(unique_transcripts)} missing."
        )
        # remove the missing transcripts by only keeping the data for the transcripts present in the mapping
        transcript_ids_intersect = list(set(unique_transcripts).intersection(set(transcript_gene_map["transcript_id"])))
        transcript_ids_intersect_boolean = np.isin(transcript_ids, transcript_ids_intersect)
        transcript_data = transcript_data.isel(
            transcript_id=transcript_ids_intersect_boolean,
            drop=True,
        )
        transcript_ids = transcript_data.coords["transcript_id"].values
        transcript_gene_map = transcript_gene_map[transcript_gene_map["transcript_id"].isin(transcript_ids_intersect)]

    # add the corresponding gene to the transcript-level expression
    log(25, "Matching gene_ids.")
    transcript_gene_dict = transcript_gene_map.set_index("transcript_id")["gene_id"].to_dict()
    gene_ids_raw = transcript_data["transcript_id"].to_series().map(transcript_gene_dict).values
    gene_ids = np.repeat(gene_ids_raw, transcript_data["abundance"].shape[1])

    # remove the transcript_id coordinate
    transcript_data = transcript_data.drop_vars("transcript_id")
    transcript_data = transcript_data.assign_coords(gene_id=gene_ids_raw)

    # rename the first dimension to gene
    transcript_data = transcript_data.rename({"transcript_id": "gene_id"})

    # get the unique genes but keep the order
    unique_genes = list(pd.Series(gene_ids).unique())

    log(25, "Creating gene abundance.")
    abundance_gene = xr.DataArray(
        transcript_data.groupby("gene_id").sum()["abundance"],
        dims=["gene_id", "file"],
    )

    log(25, "Creating gene counts.")
    counts_gene = xr.DataArray(
        transcript_data.groupby("gene_id").sum()["counts"],
        dims=["gene_id", "file"],
    )

    log(25, "Creating lengths.")
    transcript_data["abundance_length_product"] = transcript_data["abundance"] * transcript_data["length"]
    abundance_weighted_length = transcript_data.groupby("gene_id").sum()["abundance_length_product"]
    length = xr.DataArray(abundance_weighted_length / abundance_gene.data, dims=["gene_id", "file"], name="length")

    log(25, "Replacing missing lengths.")
    average_transcript_length_across_samples = transcript_data["length"].mean(axis=1)
    average_gene_length = average_transcript_length_across_samples.groupby("gene_id").mean()
    length = replace_missing_average_transcript_length(length, average_gene_length)

    # convert the counts to the desired count type
    if counts_from_abundance is not None:
        log(25, "Recreating gene counts from abundances.")
        counts_gene = convert_abundance_to_counts(
            counts_gene,
            abundance_gene,
            length,
            counts_from_abundance,
        )

    # convert to gene-level expression
    log(25, "Creating gene expression dataset.")
    gene_expression = xr.Dataset(
        data_vars={
            "abundance": abundance_gene,
            "counts": counts_gene,
            "length": length,
        },
        coords={
            "gene_id": unique_genes,
            "file_path": transcript_data.coords["file_path"].values,
        },
    )

    return gene_expression
