from logging import warning
from pathlib import Path
from typing import Union

import pandas as pd

from ..definitions import TranscriptData
from ..utils._convert_counts_to_tpm import convert_counts_to_tpm


def read_tsv(
    file_path: Union[str, Path],
    id_column: str,
    counts_column: str,
    length_column: str,
    abundance_column: str,
) -> TranscriptData:
    """Read a quantification file in tsv format.

    Args:
        file_path (Union[str, Path]): The path to the quantification file.

    Returns:
        TranscriptData: The transcript-level expression.
    """
    if not isinstance(file_path, Path):
        file_path = Path(file_path)

    if not file_path.exists():
        raise ImportError(f"The file does not exist: {file_path}")

    # read the quantification file as a tsv, tab separated and the first line is the column names
    transcript_data = pd.read_table(file_path, header=0)

    # check that the columns are in the table
    assert id_column in transcript_data.columns, f"Could not find the transcript id column `{id_column}`."
    assert counts_column in transcript_data.columns, f"Could not find the counts column `{counts_column}`."
    assert length_column in transcript_data.columns, f"Could not find the length column `{length_column}`."

    # calculate the transcript-level TPM if the abundance was not included
    if abundance_column is None:
        warning("Abundance column not provided, calculating TPM.", UserWarning)
        abundance = convert_counts_to_tpm(
            counts=transcript_data[counts_column].astype("float64").values,
            length=transcript_data[length_column].astype("float64").values,
        )
    else:
        assert abundance_column in transcript_data.columns, f"Could not find the abundance column `{abundance_column}`."
        abundance = transcript_data[abundance_column].astype("float64").values

    # create a DataFrame with the transcript-level expression
    transcripts = TranscriptData(
        transcript_id=transcript_data[id_column].values,
        counts=transcript_data[counts_column].astype("float64").values,
        length=transcript_data[length_column].astype("float64").values,
        abundance=abundance,
        inferential_replicates=None,
    )

    # return the transcript-level expression
    return transcripts
