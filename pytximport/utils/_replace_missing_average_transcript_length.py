import numpy as np
import xarray as xr


def replace_missing_average_transcript_length(
    length: xr.DataArray,
    length_gene_mean: xr.DataArray,
) -> xr.DataArray:
    """Replace missing mean transcript length at the sample level with the gene mean across samples.

    Args:
        length (xr.DataArray): The average length of transcripts at the gene level with a sample dimension.
        length_gene_mean (xr.DataArray): The mean length of the transcripts of the genes across samples.

    Returns:
        xr.DataArray: The average length of transcripts at the gene level with a sample dimension.
    """
    # Find the locations where values are missing and identify the corresponding rows
    is_nan = length.isnull()
    nan_rows = is_nan.any(dim="file")

    # Copy the length array to fill it
    length_filled = length.copy()

    # Fill genes where all values are NaN with the mean gene length
    all_nan_mask = is_nan.all(dim="file")
    length_filled = xr.where(all_nan_mask, length_gene_mean, length_filled)

    # Identify rows and lenghts with partial NaNs
    genes_with_partial_nan = length["gene_id"].where(nan_rows & ~all_nan_mask, drop=True)
    partial_nan_lengths = length.sel(gene_id=genes_with_partial_nan)

    # Calculate the geometric mean of the non-NaN values for each gene with partial NaNs
    geometric_means = xr.DataArray(
        np.exp(np.nanmean(np.log(partial_nan_lengths), axis=1)),
        dims=["gene_id"],
        coords={"gene_id": genes_with_partial_nan},
    )

    # Update `length_filled` with the partial NaN-filled values
    length_filled.loc[{"gene_id": genes_with_partial_nan}] = xr.where(
        is_nan.sel(gene_id=genes_with_partial_nan),
        geometric_means.sel(gene_id=partial_nan_lengths.gene_id),
        partial_nan_lengths,
    )

    return length_filled
