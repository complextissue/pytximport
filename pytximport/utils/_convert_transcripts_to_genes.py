from logging import log, warning
from typing import List, Literal, Optional, Union

import numpy as np
import pandas as pd
import xarray as xr

from ._convert_abundance_to_counts import convert_abundance_to_counts
from ._remove_transcript_version import remove_transcript_version
from ._replace_missing_average_transcript_length import (
    replace_missing_average_transcript_length,
)


def convert_transcripts_to_genes(
    transcript_data: xr.Dataset,
    transcript_gene_map: pd.DataFrame,
    counts_from_abundance: Optional[Literal["scaled_tpm", "length_scaled_tpm"]] = None,
) -> xr.Dataset:
    """Convert transcript-level expression to gene-level expression.

    Args:
        transcript_data (xr.Dataset): The transcript-level expression data from multiple samples.
        transcript_gene_map (pd.DataFrame): The mapping from transcripts to genes. Contains two columns: `transcript_id`
            and `gene_id`.
        counts_from_abundance (Optional[Literal["scaled_tpm", "length_scaled_tpm"]], optional): The type of counts to
            convert to. Defaults to "length_scaled_tpm".

    Returns:
        xr.Dataset: The gene-level expression data from multiple samples.
    """
    transcript_ids: Union[np.ndarray, List[str]] = transcript_data.coords["transcript_id"].values
    unique_transcripts = list(set(transcript_ids))

    # Avoid duplicates in the mapping
    transcripts_duplicated = transcript_gene_map["transcript_id"].duplicated()
    assert not any(transcripts_duplicated), "The mapping contains duplicates."

    # Check that at least one transcript is present in the mapping
    assert any(transcript_gene_map["transcript_id"].isin(unique_transcripts)), (
        "No transcripts are present in the mapping. "
        "Please make sure you are using the correct mapping and that the transcript IDs match the mapping. "
        "You may want to remove bars or transcript versions from the transcript IDs."
    )

    # Check whether there are any missing transcripts, and if so, warn the user and remove them
    if not set(unique_transcripts).issubset(set(transcript_gene_map["transcript_id"])):
        warning(
            "Not all transcripts are present in the mapping."
            + f" {len(set(unique_transcripts) - set(transcript_gene_map['transcript_id']))}"
            + f" out of {len(unique_transcripts)} missing. Removing the missing transcripts."
        )
        # Remove the missing transcripts by only keeping the data for the transcripts present in the mapping
        transcript_ids_intersect = list(set(unique_transcripts).intersection(set(transcript_gene_map["transcript_id"])))
        transcript_data = transcript_data.isel(
            transcript_id=np.isin(transcript_ids, transcript_ids_intersect),
            drop=True,
        )
        # transcript_ids = transcript_data.coords["transcript_id"].values
        transcript_gene_map = transcript_gene_map[transcript_gene_map["transcript_id"].isin(transcript_ids_intersect)]

    # Add the corresponding gene to the transcript-level expression
    log(25, "Matching gene_ids.")
    gene_ids_raw = (
        transcript_data["transcript_id"]
        .to_series()
        .map(transcript_gene_map.set_index("transcript_id")["gene_id"])
        .values
    )

    # Remove the transcript_id coordinate and rename the variable to gene_id
    transcript_data = (
        transcript_data.drop_vars("transcript_id")
        .assign_coords(gene_id=gene_ids_raw)
        .rename({"transcript_id": "gene_id"})
    )

    # Get the unique genes but keep the order
    unique_genes = pd.Series(gene_ids_raw).unique()

    log(25, "Creating gene abundance.")
    # We already calculate the abundance length product here so that we can reuse the sum
    transcript_data["abundance_length_product"] = xr.apply_ufunc(
        np.multiply,
        transcript_data["abundance"],
        transcript_data["length"],
    )
    transcript_data_summed_by_gene = transcript_data.groupby("gene_id").sum()
    abundance_gene = xr.DataArray(
        transcript_data_summed_by_gene["abundance"],
        dims=["gene_id", "file"],
    )

    log(25, "Creating gene counts.")
    counts_gene = xr.DataArray(
        transcript_data_summed_by_gene["counts"],
        dims=["gene_id", "file"],
    )

    inferential_replicates_gene = None
    if "inferential_replicates" in transcript_data.data_vars:
        log(25, "Creating inferential replicates.")
        inferential_replicates_gene = xr.DataArray(
            transcript_data_summed_by_gene["inferential_replicates"],
            dims=["gene_id", "bootstraps", "file"],
        )

    variances_gene = None
    if "variance" in transcript_data.data_vars and inferential_replicates_gene is not None:
        log(25, "Creating variances.")
        variances_gene = inferential_replicates_gene.var(dim="bootstraps", ddof=1)

    log(25, "Creating lengths.")
    length = xr.DataArray(
        transcript_data_summed_by_gene["abundance_length_product"] / abundance_gene.data,
        dims=["gene_id", "file"],
        name="length",
    )

    log(25, "Replacing missing lengths.")
    length = replace_missing_average_transcript_length(
        length,
        # Average gene length across samples
        transcript_data["length"].mean(axis=1).groupby("gene_id").mean(),
    )

    # Convert the counts to the desired count type
    if counts_from_abundance is not None:
        log(25, "Recreating gene counts from abundances.")
        counts_gene = convert_abundance_to_counts(
            counts_gene,
            abundance_gene,
            length,
            counts_from_abundance,
        )

    # Convert to gene-level expression
    log(25, "Creating gene expression dataset.")
    data_vars = {
        "abundance": abundance_gene,
        "counts": counts_gene,
        "length": length,
    }

    if inferential_replicates_gene is not None:
        data_vars["inferential_replicates"] = inferential_replicates_gene

        if variances_gene is not None:
            data_vars["variance"] = variances_gene

    gene_expression = xr.Dataset(
        data_vars=data_vars,
        coords={
            "gene_id": unique_genes,
            "file_path": transcript_data.coords["file_path"].values,
        },
    )

    return gene_expression
