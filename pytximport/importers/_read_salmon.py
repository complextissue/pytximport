from pathlib import Path
from typing import Union

import pandas as pd

from ..definitions import TranscriptData
from ._read_tsv import read_tsv


def read_salmon(
    file_path: Union[str, Path],
    id_column: str = "Name",
    counts_column: str = "NumReads",
    length_column: str = "EffectiveLength",
    abundance_column: str = "TPM",
) -> TranscriptData:
    """Read a salmon quantification file.

    Args:
        file_path (Union[str, Path]): The path to the quantification file.

    Returns:
        TranscriptData: The transcript-level expression.
    """
    if not isinstance(file_path, Path):
        file_path = Path(file_path)

    if file_path.is_dir():
        # add quant.sf to the file path
        file_path = file_path / "quant.sf"

    # check that we are importing a .sf file
    if not file_path.suffix == ".sf":
        raise ImportError("Only .sf files are supported.")

    return read_tsv(
        file_path,
        id_column=id_column,
        counts_column=counts_column,
        length_column=length_column,
        abundance_column=abundance_column,
    )
