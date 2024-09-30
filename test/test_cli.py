"""Test the cli function."""

import os
from pathlib import Path
from time import time
from typing import List


def test_cli(
    salmon_multiple_files: List[Path],
    transcript_gene_mapping_path_mouse: Path,
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
