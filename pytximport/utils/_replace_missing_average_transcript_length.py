import numpy as np
import xarray as xr


def replace_missing_average_transcript_length(
    length: xr.DataArray,
    length_gene_mean: xr.DataArray,
) -> xr.DataArray:
    """Replace missing mean transcript length at the sample level with the gene mean across samples.

    Args:
        length (xr.DataArray): The average length of transcripts at the gene level with a sample dimension.
        length_gene_mean (xr.DataArray): The mean length of the transcript of the genes across samples.

    Returns:
        xr.DataArray: The average length of transcripts at the gene level with a sample dimension.
    """
    # get the rows of the DataArray with missing values
    nan_rows = np.where(length.isnull().any(dim="file") == True)[0]  # noqa: E712

    for nan_idx in nan_rows:
        row = length.isel({"gene_id": nan_idx})
        gene_id = row.coords["gene_id"].data
        row_data = row.data
        column_is_nan = np.isnan(row_data)

        if len(row_data) == np.sum(column_is_nan):
            # get the average length of the gene
            average_gene_length = length_gene_mean.loc[{"gene_id": gene_id}].data

            # if the gene is not present in the gene mean, replace with 0
            if np.isnan(average_gene_length):
                average_gene_length = 0

        else:
            # replace with the geometric mean
            average_gene_length = np.exp(np.mean(np.log(length.loc[{"gene_id": gene_id}].data[~column_is_nan])))

        # replace the missing row with the average gene length
        length.loc[{"gene_id": gene_id}] = length.loc[{"gene_id": gene_id}].fillna(average_gene_length)

    return length
