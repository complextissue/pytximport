import numpy as np
from numpy.typing import NDArray


def convert_counts_to_tpm(
    counts: NDArray,
    length: NDArray,
) -> NDArray:
    """Convert transcript-level counts to TPM.

    Args:
        counts (NDArray): The transcript-level counts.
        length (NDArray): The length of the transcripts.

    Returns:
        NDArray: The transcript-level expression data with the TPM.
    """
    return np.array(counts * (1e6 / np.sum(counts / length)) / length)
