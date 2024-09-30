import importlib.util
from logging import warning
from pathlib import Path
from typing import Literal, Optional, Union

import numpy as np
import pandas as pd
from numpy.typing import ArrayLike

from ..definitions import TranscriptData
from ..utils._convert_counts_to_tpm import convert_counts_to_tpm


def parse_dataframe(
    transcript_dataframe: pd.DataFrame,
    id_column: str,
    counts_column: str,
    length_column: str,
    abundance_column: Optional[str] = None,
    recompute_counts: bool = False,
) -> TranscriptData:
    """Parse a DataFrame with the transcript-level expression.

    Args:
        transcript_dataframe (pd.DataFrame): The DataFrame with the transcript-level expression.
        id_column (str): The column name for the transcript id.
        counts_column (str): The column name for the counts.
        length_column (str): The column name for the length.
        abundance_column (Optional[str], optional): The column name for the abundance. Defaults to None.
        recompute_counts (bool, optional): Whether inferential replicates will be used to recompute counts and
            abundances. If true, the counts and abundances will not be read from the file. Defaults to False.

    Returns:
        TranscriptData: The transcript-level expression.
    """
    # Check that the columns are in the table
    assert id_column in transcript_dataframe.columns, f"Could not find the transcript id column `{id_column}`."
    assert length_column in transcript_dataframe.columns, f"Could not find the length column `{length_column}`."

    counts: Optional[ArrayLike]
    if not recompute_counts:
        assert counts_column in transcript_dataframe.columns, f"Could not find the counts column `{counts_column}`."

        # Calculate the transcript-level TPM if the abundance was not included
        if abundance_column is None:
            warning("Abundance column not provided, calculating TPM.", UserWarning)
            abundance = convert_counts_to_tpm(
                counts=transcript_dataframe[counts_column].values,  # type: ignore
                length=transcript_dataframe[length_column].values,  # type: ignore
            )
        else:
            assert (
                abundance_column in transcript_dataframe.columns
            ), f"Could not find the abundance column `{abundance_column}`."
            abundance = transcript_dataframe[abundance_column].values  # type: ignore

        counts = transcript_dataframe[counts_column].values  # type: ignore
    else:
        counts = None
        abundance = None

    # Create a DataFrame with the transcript-level expression
    transcripts = TranscriptData(
        transcript_id=transcript_dataframe[id_column].values,  # type: ignore
        counts=counts,
        length=transcript_dataframe[length_column].values,  # type: ignore
        abundance=abundance,
        inferential_replicates=None,
    )

    # Return the transcript-level expression
    return transcripts


def read_tsv(
    file_path: Union[str, Path],
    id_column: str,
    counts_column: str,
    length_column: str,
    abundance_column: Optional[str] = None,
    recompute_counts: bool = False,
) -> TranscriptData:
    """Read a quantification file in tsv format.

    Args:
        file_path (Union[str, Path]): The path to the quantification file.
        id_column (str): The column name for the transcript id.
        counts_column (str): The column name for the counts.
        length_column (str): The column name for the length.
        abundance_column (Optional[str], optional): The column name for the abundance. Defaults to None.
        recompute_counts (bool, optional): Whether inferential replicates will be used to recompute counts and
            abundances. If true, the counts and abundances will not be read from the file. Defaults to False.

    Returns:
        TranscriptData: The transcript-level expression.
    """
    if not isinstance(file_path, Path):
        file_path = Path(file_path)

    if not file_path.exists():
        raise ImportError(f"The file does not exist: {file_path}")

    # Read the quantification file as a tsv, tab separated with the first line being the column names
    usecols = [id_column, length_column]
    dtype = {id_column: str, length_column: np.float64}

    if not recompute_counts:
        usecols.append(counts_column)
        dtype[counts_column] = np.float64

        if abundance_column is not None:
            usecols.append(abundance_column)
            dtype[abundance_column] = np.float64

    # Check if pyarrow is available
    engine: Literal["pyarrow", "c"] = "pyarrow" if importlib.util.find_spec("pyarrow") is not None else "c"

    if engine != "pyarrow":
        warning("pyarrow is not available, consider installing it to improve import performance.")

    transcript_dataframe = pd.read_table(
        file_path,
        header=0,
        sep="\t",
        compression=("gzip" if file_path.suffix == ".gz" else None),
        engine=engine,
        usecols=usecols,
        dtype=dtype,
        na_filter=False,
    )

    return parse_dataframe(
        transcript_dataframe,
        id_column=id_column,
        counts_column=counts_column,
        length_column=length_column,
        abundance_column=abundance_column,
        recompute_counts=recompute_counts,
    )
