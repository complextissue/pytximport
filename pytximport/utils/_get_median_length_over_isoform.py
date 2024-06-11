import numpy as np
import pandas as pd
import xarray as xr

from ._remove_transcript_version import remove_transcript_version


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

    # get the transcript lengths
    transcript_gene_dict = transcript_gene_map.set_index("transcript_id")["gene_id"].to_dict()
    gene_ids = transcript_data["transcript_id"].to_series().map(transcript_gene_dict).values

    # check that the gene ids are not all nan
    assert not any(pd.isna(gene_ids)), "Not all transcript ids could be mapped to gene ids. Please check the mapping."

    transcript_data_copy = transcript_data.drop_vars("transcript_id")
    transcript_data_copy = transcript_data_copy.assign_coords(gene_id=gene_ids)
    transcript_data_copy = transcript_data_copy.rename({"transcript_id": "gene_id"})

    # get the row mean across samples for each transcript
    average_transcript_length_across_samples = transcript_data_copy["length"].mean(axis=1)
    median_gene_length = average_transcript_length_across_samples.groupby("gene_id").median().to_dataframe()

    transcript_median_gene_length = [median_gene_length.loc[gene_id, "length"] for gene_id in gene_ids]
    transcript_median_gene_length_repeated = np.reshape(
        np.repeat(
            transcript_median_gene_length,
            transcript_data["abundance"].shape[1],
        ),
        transcript_data["abundance"].shape,
    )

    transcript_data["median_isoform_length"] = xr.DataArray(
        transcript_median_gene_length_repeated,
        dims=("transcript_id", "file"),
    )

    return transcript_data
