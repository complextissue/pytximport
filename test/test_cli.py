"""Test the cli function."""

import os
from pathlib import Path
from time import time
from typing import List

import pandas as pd


def test_cli(
    salmon_multiple_files: List[Path],
    transcript_gene_mapping_path_mouse: Path,
    annotation_file: Path,
) -> None:
    """Test the cli function."""
    file_paths = " -i ".join([str(file) for file in salmon_multiple_files])

    # Timestamp
    current_time = int(time())

    # Check whether the current directory can be written to
    if not os.access(os.getcwd(), os.W_OK):
        raise PermissionError("The current directory cannot be written to.")

    command = (
        f"pytximport -i {file_paths} -t salmon -m {str(transcript_gene_mapping_path_mouse)}"
        f" -o ./pytximport_cli_test_{current_time}.csv"
    )

    os.system(command)  # nosec

    # Check that the output file was created
    assert os.path.exists(f"./pytximport_cli_test_{current_time}.csv"), "Output file was not created."

    # Remove the temporary file
    os.system(f"rm ./pytximport_cli_test_{current_time}.csv")  # nosec

    # Test saving as h5ad
    command_ad = f"{command.replace('.csv', '.h5ad')} --output_format h5ad"

    os.system(command_ad)  # nosec

    # Check that the output file was created
    assert os.path.exists(f"./pytximport_cli_test_{current_time}.h5ad"), "AnnData output file was not created."

    # Remove the temporary file
    os.system(f"rm ./pytximport_cli_test_{current_time}.h5ad")  # nosec

    # Test creating a transcript-to-gene mapping
    command = (
        f"pytximport create-map -i {str(annotation_file)} -o ./pytximport_cli_test_map_{current_time}.csv -ow"
        " --source-field transcript_id --target-field gene_id --target-field gene_biotype"
    )

    os.system(command)  # nosec

    # Check that the output file was created
    assert os.path.exists(f"./pytximport_cli_test_map_{current_time}.csv"), "Mapping file was not created."

    transcript_gene_map = pd.read_csv(f"./pytximport_cli_test_map_{current_time}.csv", header=0, index_col=None)

    # Check that the columns are present
    assert "transcript_id" in transcript_gene_map.columns, "The transcript_id column is not present."
    assert "gene_id" in transcript_gene_map.columns, "The gene_id column is not present."
    assert "gene_biotype" in transcript_gene_map.columns, "The gene_biotype column is not present."

    # Remove the temporary file
    os.system(f"rm ./pytximport_cli_test_map_{current_time}.csv")  # nosec
