import numpy as np
import pandas as pd
import xarray as xr


def get_median_length_over_isoform(
    transcript_data: xr.Dataset,
    transcript_gene_map: pd.DataFrame,
) -> xr.Dataset:
    """Get the median length of the gene over all isoforms.

    Args:
        length (xr.Dataset): The transcript data containing the length of the transcripts.
        transcript_gene_map (pd.DataFrame): The mapping of transcripts to genes.
        ignore_after_bar (bool, optional): Whether to ignore the part of the transcript ID after the bar.
            Defaults to True.

    Returns:
        xr.Dataset: The updated transcript data with the median gene length contained in the `median_isoform_length`
            variable.
    """
    assert "length" in transcript_data.data_vars, "The transcript data does not contain a `length` variable."

    # Get the gene ids for each transcript
    gene_ids = (
        transcript_data["transcript_id"]
        .to_series()
        .map(transcript_gene_map.set_index("transcript_id")["gene_id"].to_dict())
        .values
    )

    # Check that no gene ids is nan
    assert not any(pd.isna(gene_ids)), "Not all transcript ids could be mapped to gene ids. Please check the mapping."

    # Get the row mean across samples for each transcript
    median_gene_length = (
        transcript_data.drop_vars("transcript_id")
        .assign_coords(gene_id=gene_ids)
        .rename({"transcript_id": "gene_id"})["length"]
        .mean(dim="file")
        .groupby("gene_id")
        .median()
        .to_dataframe()
    )

    transcript_data["median_isoform_length"] = xr.DataArray(
        np.reshape(
            np.repeat(
                pd.Series(gene_ids).map(median_gene_length["length"]).to_numpy(),
                transcript_data["abundance"].shape[1],
            ),
            transcript_data["abundance"].shape,
        ),
        dims=("transcript_id", "file"),
    )

    return transcript_data
