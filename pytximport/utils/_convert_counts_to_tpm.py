import numpy as np


def convert_counts_to_tpm(
    counts: np.ndarray,
    length: np.ndarray,
) -> np.ndarray:
    """Convert transcript-level counts to TPM.

    Args:
        counts (np.ndarray): The transcript-level counts.
        length (np.ndarray): The length of the transcripts.

    Returns:
        np.ndarray: The transcript-level expression data with the TPM.
    """
    normalization_factor = 1e6 / np.sum(counts / length)
    return np.array(counts * normalization_factor / length)
