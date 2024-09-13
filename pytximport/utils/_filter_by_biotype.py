from logging import log
from typing import List

import numpy as np
import xarray as xr


def filter_by_biotype(
    transcript_data: xr.Dataset,
    biotype_filter: List[str],
    id_column: str = "transcript_id",
) -> xr.Dataset:
    """Filter the transcripts by biotype.

    This function filters the transcripts by biotype. The biotype is assumed to be present in the transcript_id
    separated by a bar. The biotype is checked against the biotype_filter and the transcripts that match the biotype
    are kept. This function is provided mainly for internal use if `biotype_filter` is provided to the main function.

    Args:
        transcript_data (xr.Dataset): The transcript-level expression data from multiple samples.
        biotype_filter (List[str]): The biotypes to filter the transcripts by.
        id_column (str, optional): The column name for the transcript ID. Defaults to "transcript_id".

    Returns:
        xr.Dataset: The transcript-level expression data from multiple samples with the transcripts filtered by biotype.
    """
    transcript_ids = transcript_data.coords[id_column].values

    # This only works if the transcript_id contains the biotype as one of the bar-separated fields
    assert any(
        "|" in transcript_id for transcript_id in transcript_ids
    ), "The transcript_id does not contain the biotype."

    transcript_id_fields = [transcript_id.split("|") for transcript_id in transcript_ids]

    # Iterate through the biotypes and mark a transcript to be kept if a biotype is a match
    transcript_keep_boolean = np.zeros(len(transcript_id_fields))

    for biotype in biotype_filter:
        transcript_keep_boolean = np.logical_or(
            transcript_keep_boolean,
            [(biotype in transcript_id_field) for transcript_id_field in transcript_id_fields],
        )

    # Check that at least one transcript is protein-coding
    assert any(transcript_keep_boolean), "No transcripts with the desired biotype are present in the data."

    # Calculate the total abundance before filtering
    total_abundance = transcript_data["abundance"].sum(axis=0)

    transcript_data = transcript_data.isel(
        transcript_id=transcript_keep_boolean,
        drop=True,
    )
    log(
        25,
        f"Removed {len(transcript_keep_boolean) - sum(transcript_keep_boolean)} transcripts with other biotypes.",
    )
    transcript_ids = transcript_data.coords[id_column].values

    # Recalculate the abundance for each sample
    log(25, "Recalculating the abundance after filtering by biotype.")
    new_abundance = transcript_data["abundance"].sum(axis=0)
    ratio = total_abundance / new_abundance
    transcript_data["abundance"] = (transcript_data["abundance"].T * ratio).T

    return transcript_data
