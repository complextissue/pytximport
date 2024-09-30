from logging import warning
from pathlib import Path
from typing import Optional, Union

import numpy as np
from h5py import File

from ..definitions import InferentialReplicates, TranscriptData
from ..utils._convert_counts_to_tpm import convert_counts_to_tpm
from ._read_tsv import read_tsv


def read_inferential_replicates_kallisto(
    file_path: Union[str, Path],
) -> InferentialReplicates:
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
        raise ImportError("Only .h5 files are supported.")

    if not file_path.exists():
        raise ImportError(f"The file does not exist: {file_path}")

    with File(file_path, "r") as f:
        if "/bootstrap" not in f.keys():
            raise ImportError("No inferential replicates found.")

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


def read_kallisto(
    file_path: Union[str, Path],
    id_column: str = "aux/ids",
    counts_column: str = "est_counts",
    length_column: str = "aux/eff_lengths",
    abundance_column: Optional[str] = None,
    inferential_replicates: bool = False,
    recompute_counts: bool = False,
) -> TranscriptData:
    """Read a kallisto quantification file.

    Args:
        file_path (Union[str, Path]): The path to the quantification file.
        id_column (str, optional): The column name for the transcript id. Defaults to "aux/ids".
        counts_column (str, optional): The column name for the counts. Defaults to "est_counts".
        length_column (str, optional): The column name for the length. Defaults to "aux/eff_lengths".
        abundance_column (Optional[str], optional): The column name for the abundance. Defaults to None.
        inferential_replicates (bool, optional): Whether to read inferential replicates. Defaults to False.
        recompute_counts (bool, optional): Whether inferential replicates will be used to recompute counts and
            abundances. If true, the counts and abundances will not be read from the file. Defaults to False.

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

    # Read the quantification file
    if file_path.suffix == ".h5":
        with File(file_path, "r") as f:
            # Get the transcript-level expression
            transcript_ids = f.file[id_column][:].astype(str)
            counts = f.file[counts_column][:]
            length = f.file[length_column][:]

            # Get the abundance if it was specified, else it will be calculated
            if abundance_column is not None:
                abundance = f.file[abundance_column][:]

        # Check that the length of the counts, length, and abundances are the same
        assert (
            len(transcript_ids) == len(counts) == len(length)
        ), "The transcript ids, counts and length have different length."

        # Calculate the transcript-level TPM if the abundance was not included
        if abundance_column is None:
            warning("Abundance column not provided, calculating TPM.")
            abundance = convert_counts_to_tpm(counts, length)
        else:
            assert len(transcript_ids) == len(abundance), "The transcript ids and abundance have different length."

        # Create a DataFrame with the transcript-level expression
        transcripts = TranscriptData(
            transcript_id=transcript_ids,
            counts=counts,
            length=length,
            abundance=abundance,
            inferential_replicates=None,
        )

        if inferential_replicates:
            transcripts["inferential_replicates"] = read_inferential_replicates_kallisto(
                file_path,
            )

    elif file_path.suffix == ".tsv":
        transcripts = read_tsv(
            file_path,
            id_column=id_column,
            counts_column=counts_column,
            length_column=length_column,
            abundance_column=abundance_column,
            recompute_counts=recompute_counts,
        )

        if inferential_replicates:
            # Check whether there is a matching inferential replicates file
            inferential_replicates_file = file_path.parent / "abundance.h5"

            if inferential_replicates_file.exists():
                transcripts["inferential_replicates"] = read_inferential_replicates_kallisto(
                    inferential_replicates_file,
                )
            else:
                raise ImportError("No inferential replicates found. Please provide the path to the abundance.h5 file.")

    return transcripts
