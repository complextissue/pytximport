from logging import warning
from pathlib import Path
from typing import Union

from ..definitions import TranscriptData
from ._read_tsv import read_tsv


def read_rsem(
    file_path: Union[str, Path],
    id_column: str = "transcript_id",  # Or "gene_id" for gene-level quantification
    counts_column: str = "expected_count",
    length_column: str = "effective_length",
    abundance_column: str = "TPM",
    gene_level: bool = False,
) -> TranscriptData:
    """Read an RSEM quantification file.

    Args:
        file_path (Union[str, Path]): The path to the quantification file.
        id_column (str, optional): The column name for the transcript ID. Defaults to "transcript_id".
        counts_column (str, optional): The column name for the counts. Defaults to "expected_count".
        length_column (str, optional): The column name for the length. Defaults to "effective_length".
        abundance_column (str, optional): The column name for the abundance. Defaults to "TPM".
        gene_level (bool, optional): Whether the quantification is at the gene level. Defaults to False.

    Returns:
        TranscriptData: The transcript-level expression.
    """
    if not isinstance(file_path, Path):
        file_path = Path(file_path)

    if file_path.is_dir():
        if gene_level:
            file_paths = list(file_path.glob("*.genes.results.gz"))
            file_identifier = "genes.results.gz"
        else:
            file_paths = list(file_path.glob("*.isoforms.results.gz"))
            file_identifier = "isoforms.results.gz"

        if not file_paths:
            raise ImportError(f"No {file_identifier} files found in the directory.")

        if len(file_paths) > 1:
            raise ImportError(f"Multiple {file_identifier} files found in the directory.")

        file_path = file_paths[0]

    # Check that we are importing a .sf file
    if not file_path.suffix == ".gz" and not file_path.suffix == ".results":
        raise ImportError("Only .gz and .results files are supported.")

    if gene_level and "transcript" in id_column:
        warning("Gene-level quantification file with transcript-level column name. Please check the column name.")
    elif not gene_level and "gene" in id_column:
        warning("Transcript-level quantification file with gene-level column name. Please check the column name.")

    transcript_data = read_tsv(
        file_path,
        id_column=id_column,
        counts_column=counts_column,
        length_column=length_column,
        abundance_column=abundance_column,
    )

    # Set the minimum length to 1
    transcript_data["length"] = transcript_data["length"].clip(min=1)  # type: ignore

    return transcript_data
