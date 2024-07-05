from pathlib import Path
from typing import Literal, Optional, Union
from warnings import warn

import numpy as np
import pandas as pd
from h5py import File

from ..definitions import InferentialReplicates, TranscriptData
from ..utils._convert_counts_to_tpm import convert_counts_to_tpm


def read_inferential_replicates_kallisto(
    file_path: Union[str, Path],
) -> Union[InferentialReplicates, None]:
    """Read inferential replicates from a kallisto quantification file.

    Args:
        file_path (Union[str, Path]): The path to the quantification file.

    Returns:
        InferentialReplicates: The inferential replicates.
    """
    if not isinstance(file_path, Path):
        file_path = Path(file_path)

    if file_path.is_dir():
        file_path = file_path / "abundance.h5"

    if file_path.suffix != ".h5":
        return None

    if not file_path.exists():
        return None

    with File(file_path, "r") as f:
        if "/bootstrap" in f.keys():
            bootstraps = f["/bootstrap"]
            target_count = len(bootstraps["bs0"])
            inferential_replicates = np.zeros((target_count, len(bootstraps.keys())))

            for bootstrap_idx in range(len(bootstraps.keys())):
                inferential_replicates[:, bootstrap_idx] = bootstraps[f"bs{bootstrap_idx}"][:]

            vars = np.var(inferential_replicates, axis=1, ddof=1)
            return InferentialReplicates(
                variance=vars,
                replicates=inferential_replicates,
            )
        else:
            return None


def read_kallisto(
    file_path: Union[str, Path],
    id_column: str = "aux/ids",
    counts_column: str = "est_counts",
    length_column: str = "aux/eff_lengths",
    abundance_column: Optional[str] = None,
) -> TranscriptData:
    """Read a kallisto quantification file.

    Args:
        file_path (Union[str, Path]): The path to the quantification file.

    Returns:
        TranscriptData: The transcript-level expression.
    """
    if not isinstance(file_path, Path):
        file_path = Path(file_path)

    if file_path.is_dir():
        file_path = file_path / "abundance.h5"

    if file_path.suffix not in [".h5", ".tsv"]:
        raise ImportError("Only .h5 and .tsv files are supported.")

    if not file_path.exists():
        raise ImportError(f"The file does not exist: {file_path}")

    # read the quantification file
    if file_path.suffix == ".h5":
        with File(file_path, "r") as f:
            # get the transcript-level expression
            transcript_ids = f.file[id_column][:].astype(str)
            counts = f.file[counts_column][:]
            length = f.file[length_column][:]

            # get the abundance if it was specified, else it will be calculated
            if abundance_column is not None:
                abundance = f.file[abundance_column][:]

    elif file_path.suffix == ".tsv":
        # read the quantification file as a tsv, tab separated and the first line is the column names
        transcript_data = pd.read_table(file_path, header=0)

        # check that the columns are in the table
        assert id_column in transcript_data.columns, f"Could not find the transcript id column `{id_column}`."
        assert counts_column in transcript_data.columns, f"Could not find the counts column `{counts_column}`."
        assert length_column in transcript_data.columns, f"Could not find the length column `{length_column}`."

        transcript_ids = transcript_data[id_column].values
        counts = transcript_data[counts_column].astype("float64").values
        length = transcript_data[length_column].astype("float64").values

        if abundance_column is not None:
            assert (
                abundance_column in transcript_data.columns
            ), f"Could not find the abundance column `{abundance_column}`."
            abundance = transcript_data[abundance_column].astype("float64").values

    # check that the length of the counts, length, and abundances are the same
    assert (
        len(transcript_ids) == len(counts) == len(length)
    ), "The transcript ids, counts and length have different length."

    # calculate the transcript-level TPM if the abundance was not included
    if abundance_column is None:
        warn("Abundance column not provided, calculating TPM.", UserWarning)
        abundance = convert_counts_to_tpm(counts, length)
    else:
        assert len(transcript_ids) == len(abundance), "The transcript ids and abundance have different length."

    # create a DataFrame with the transcript-level expression
    transcripts = TranscriptData(
        transcript_id=transcript_ids,
        counts=counts,
        length=length,
        abundance=abundance,
        inferential_replicates=None,
    )

    transcripts["inferential_replicates"] = read_inferential_replicates_kallisto(
        file_path,
    )

    # return the transcript-level expression
    return transcripts
