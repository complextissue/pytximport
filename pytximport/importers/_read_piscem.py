from pathlib import Path
from typing import Literal, Union

import numpy as np
import pandas as pd

from ..definitions import InferentialReplicates, TranscriptData
from ._read_tsv import read_tsv


def read_inferential_replicates_piscem(
    file_path: Union[str, Path],
) -> InferentialReplicates:
    """Read inferential replicates from a piscem quantification file.

    Args:
        file_path (Union[str, Path]): The path to the quantification file. The file should be a .quant file that is
            colocated with the inferential replicates file (.infreps.pq).

    Returns:
        InferentialReplicates: The inferential replicates.
    """
    if not isinstance(file_path, Path):
        file_path = Path(file_path)

    # Add .infreps.pq to the stem
    file_path = file_path.with_suffix(".infreps.pq")

    # Check whether the file exists
    if not file_path.exists():
        raise ImportError(f"The file does not exist: {file_path}")

    # Read the inferential replicates
    try:
        # Pandas will attempt to use pyarrow by default with a fallback to fastparquet
        bootstrap_data = pd.read_parquet(file_path)
    except ImportError:
        raise ImportError(
            "Could not read inferential replicates. "
            "Either pyarrow or fastparquet is required to read inferential replicates."
        )

    return InferentialReplicates(
        variance=np.var(bootstrap_data.to_numpy(dtype=np.float64), axis=1, ddof=1),
        replicates=bootstrap_data.to_numpy(dtype=np.float64),
    )


def read_piscem(
    file_path: Union[str, Path],
    id_column: str = "target_name",
    counts_column: str = "ecount",
    length_column: str = "eeln",
    abundance_column: str = "tpm",
    inferential_replicates: bool = False,
    recompute_counts: bool = False,
) -> TranscriptData:
    """Read a piscem-infer quantification file.

    Args:
        file_path (Union[str, Path]): The path to the quantification file.
        id_column (str, optional): The column name for the transcript id. Defaults to "Name".
        counts_column (str, optional): The column name for the counts. Defaults to "NumReads".
        length_column (str, optional): The column name for the length. Defaults to "EffectiveLength".
        abundance_column (str, optional): The column name for the abundance. Defaults to "TPM".
        aux_dir_name (Literal["aux_info", "aux"], optional): The name of the aux directory. Defaults to "aux_info".
        inferential_replicates (bool, optional): Whether to read inferential replicates. Defaults to False.
        recompute_counts (bool, optional): Whether inferential replicates will be used to recompute counts and
            abundances. If true, the counts and abundances will not be read from the file. Defaults to False.

    Returns:
        TranscriptData: The transcript-level expression.
    """
    if not isinstance(file_path, Path):
        file_path = Path(file_path)

    # Check that we are importing a .quant file
    if not file_path.suffix == ".quant" and not file_path.suffix == ".gz":
        raise ImportError("Only .quant and .gz files are supported.")

    transcript_data = read_tsv(
        file_path,
        id_column=id_column,
        counts_column=counts_column,
        length_column=length_column,
        abundance_column=abundance_column,
        recompute_counts=recompute_counts,
    )

    if inferential_replicates:
        transcript_data["inferential_replicates"] = read_inferential_replicates_piscem(
            file_path,
        )

    return transcript_data
