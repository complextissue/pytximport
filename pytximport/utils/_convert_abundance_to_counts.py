from logging import log
from typing import Literal

from xarray import DataArray


def convert_abundance_to_counts(
    counts: DataArray,
    abundance: DataArray,
    length: DataArray,
    counts_from_abundance: Literal["scaled_tpm", "length_scaled_tpm"],
) -> DataArray:
    """Convert transcript-level abundance to counts, either as TPM or TPM scaled by the length.

    Args:
        counts (DataArray): The original counts.
        abundance (DataArray): The transcript-level abundance.
        length (DataArray): The length of the transcripts.
        counts_from_abundance (Literal["scaled_tpm", "length_scaled_tpm"]): The type of counts to convert to.

    Returns:
        DataArray: The transcript-level expression data with the counts.
    """
    if counts_from_abundance == "scaled_tpm":
        # set the counts to the TPM
        log(25, "Setting the counts to scaled TPM.")
        counts_transformed = abundance
    elif counts_from_abundance == "length_scaled_tpm":
        # convert the TPM to counts and scale by the length
        log(25, "Setting counts to length scaled TPM.")
        counts_transformed = abundance * length.mean(axis=1)
    else:
        raise ValueError("The count transform must be 'scaled_tpm' or 'length_scaled_tpm'.")

    # scale the counts
    column_counts = counts.sum(axis=0)
    new_counts = counts_transformed.sum(axis=0)
    ratio = column_counts / new_counts
    counts_transformed = (counts_transformed.T * ratio).T

    return counts_transformed
