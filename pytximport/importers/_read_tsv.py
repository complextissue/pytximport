from logging import warning
from pathlib import Path
from typing import Optional, Union

import pandas as pd

from ..definitions import TranscriptData
from ..utils._convert_counts_to_tpm import convert_counts_to_tpm


def parse_dataframe(
    transcript_dataframe: pd.DataFrame,
    id_column: str,
    counts_column: str,
    length_column: str,
    abundance_column: Optional[str] = None,
) -> TranscriptData:
    """Parse a DataFrame with the transcript-level expression.

    Args:
        transcript_dataframe (pd.DataFrame): The DataFrame with the transcript-level expression.
        id_column (str): The column name for the transcript id.
        counts_column (str): The column name for the counts.
        length_column (str): The column name for the length.
        abundance_column (Optional[str], optional): The column name for the abundance. Defaults to None.

    Returns:
        TranscriptData: The transcript-level expression.
    """
    # check that the columns are in the table
    assert id_column in transcript_dataframe.columns, f"Could not find the transcript id column `{id_column}`."
    assert counts_column in transcript_dataframe.columns, f"Could not find the counts column `{counts_column}`."
    assert length_column in transcript_dataframe.columns, f"Could not find the length column `{length_column}`."

    # calculate the transcript-level TPM if the abundance was not included
    if abundance_column is None:
        warning("Abundance column not provided, calculating TPM.", UserWarning)
        abundance = convert_counts_to_tpm(
            counts=transcript_dataframe[counts_column].astype("float64").values,
            length=transcript_dataframe[length_column].astype("float64").values,
        )
    else:
        assert (
            abundance_column in transcript_dataframe.columns
        ), f"Could not find the abundance column `{abundance_column}`."
        abundance = transcript_dataframe[abundance_column].astype("float64").values

    # create a DataFrame with the transcript-level expression
    transcripts = TranscriptData(
        transcript_id=transcript_dataframe[id_column].values,
        counts=transcript_dataframe[counts_column].astype("float64").values,
        length=transcript_dataframe[length_column].astype("float64").values,
        abundance=abundance,
        inferential_replicates=None,
    )

    # return the transcript-level expression
    return transcripts


def read_tsv(
    file_path: Union[str, Path],
    id_column: str,
    counts_column: str,
    length_column: str,
    abundance_column: Optional[str] = None,
) -> TranscriptData:
    """Read a quantification file in tsv format.

    Args:
        file_path (Union[str, Path]): The path to the quantification file.
        id_column (str): The column name for the transcript id.
        counts_column (str): The column name for the counts.
        length_column (str): The column name for the length.
        abundance_column (Optional[str], optional): The column name for the abundance. Defaults to None.

    Returns:
        TranscriptData: The transcript-level expression.
    """
    if not isinstance(file_path, Path):
        file_path = Path(file_path)

    if not file_path.exists():
        raise ImportError(f"The file does not exist: {file_path}")

    # read the quantification file as a tsv, tab separated and the first line is the column names
    if file_path.suffix == ".gz":
        transcript_dataframe = pd.read_table(file_path, header=0, compression="gzip", sep="\t")
    else:
        transcript_dataframe = pd.read_table(file_path, header=0, sep="\t")

    return parse_dataframe(
        transcript_dataframe,
        id_column=id_column,
        counts_column=counts_column,
        length_column=length_column,
        abundance_column=abundance_column,
    )
