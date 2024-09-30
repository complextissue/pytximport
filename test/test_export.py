"""Test the cli function."""

import os
from pathlib import Path
from time import time
from typing import List

from pytximport import tximport


def test_export(
    salmon_multiple_files: List[Path],
    transcript_gene_mapping_path_mouse: Path,
) -> None:
    """Test saving the output of the tximport function as a csv and h5ad file."""
    if not os.access(os.getcwd(), os.W_OK):
        raise PermissionError("The current directory cannot be written to.")

    for output_type in ["anndata", "xarray"]:
        current_time = int(time())

        _ = tximport(
            salmon_multiple_files,
            "salmon",
            transcript_gene_mapping_path_mouse,
            output_type=output_type,
            output_format="csv",
            output_path=f"./pytximport_cli_test_{int(time())}.csv",
        )

        # Check that the output file was created
        assert os.path.exists(f"./pytximport_cli_test_{current_time}.csv"), "Output file was not created."

        # Remove the temporary file
        os.system(f"rm ./pytximport_cli_test_{current_time}.csv")  # nosec

    # Test saving as h5ad
    _ = tximport(
        salmon_multiple_files,
        "salmon",
        transcript_gene_mapping_path_mouse,
        output_type="anndata",
        output_format="h5ad",
        output_path=f"./pytximport_cli_test_{current_time}.h5ad",
    )

    # Check that the output file was created
    assert os.path.exists(f"./pytximport_cli_test_{current_time}.h5ad"), "AnnData output file was not created."

    # Remove the temporary file
    os.system(f"rm ./pytximport_cli_test_{current_time}.h5ad")  # nosec
