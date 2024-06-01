"""Test importing salmon quantification files."""

from pathlib import Path
from typing import List

import anndata as ad
import pandas as pd
import xarray as xr

from pytximport import tximport


def test_correctness(
    fabry_disease_files: List[Path],
) -> None:
    """Test importing salmon quantification files.

    Args:
        salmon_multiple_files (Path): [description]
    """
    fabry_directory = fabry_disease_files[0].parent

    for counts_from_abundance in [None, "scaled_tpm", "length_scaled_tpm"]:
        result = tximport(
            fabry_disease_files,
            "salmon",
            fabry_directory / "transcript_gene_mapping_human.csv",
            ignore_transcript_version=True,
            ignore_after_bar=True,
            counts_from_abundance=counts_from_abundance,  # type: ignore
        )

        assert isinstance(result, xr.Dataset), "The result is not an xarray Dataset."

        df_result = pd.DataFrame(
            result["counts"].data,
            index=result.coords["gene_id"],
            columns=result.coords["file_path"],
        ).sort_index()

        # load in the comparison data generated with the R package tximport
        if counts_from_abundance is None:
            df_result_tximport = pd.read_csv(
                fabry_directory / "counts_tximport_no.csv",
                index_col=0,
                header=0,
            )
        elif counts_from_abundance == "scaled_tpm":
            df_result_tximport = pd.read_csv(
                fabry_directory / "counts_tximport_scaledTPM.csv",
                index_col=0,
                header=0,
            )
        elif counts_from_abundance == "length_scaled_tpm":
            df_result_tximport = pd.read_csv(
                fabry_directory / "counts_tximport_lengthScaledTPM.csv",
                index_col=0,
                header=0,
            )

        # since tximport does not use the file path as column, use the columns from the comparison data
        df_result.columns = df_result_tximport.columns

        # check that the data is the same
        pd.testing.assert_frame_equal(df_result, df_result_tximport)
