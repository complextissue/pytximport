"""Test importing salmon quantification files."""

from pathlib import Path
from typing import List

import numpy as np
import pandas as pd
import xarray as xr

from pytximport import tximport


def test_correctness(
    fabry_disease_files: List[Path],
) -> None:
    """Test importing salmon quantification files.

    Args:
        fabry_disease_files (List[Path]): The paths to the quantification files.
    """
    fabry_directory = fabry_disease_files[0].parent

    for counts_from_abundance in [None, "scaled_tpm", "length_scaled_tpm"]:
        result = tximport(
            fabry_disease_files,
            "salmon",
            fabry_directory / "transcript_gene_mapping_human.tsv",
            ignore_transcript_version=True,
            ignore_after_bar=True,
            output_type="xarray",
            counts_from_abundance=counts_from_abundance,  # type: ignore
        )

        assert isinstance(result, xr.Dataset), "The result is not an xarray Dataset."

        df_result = pd.DataFrame(
            result["counts"].data,
            index=result.coords["gene_id"],
            columns=result.coords["file_path"],
        ).sort_index()

        # Load in the comparison data generated with the R package tximport
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

        # Since tximport does not use the file path as column, use the columns from the comparison data
        df_result.columns = df_result_tximport.columns

        # Check that the data is the same
        pd.testing.assert_frame_equal(df_result, df_result_tximport)


def test_correctness_transcript_level(
    salmon_file: Path,
) -> None:
    """Test that the transcript level data is correct.

    Args:
        salmon_file (Path): Path to the salmon quantification file.
    """
    data_directory = salmon_file.parent.parent

    for counts_from_abundance in ["scaled_tpm", "dtu_scaled_tpm"]:
        result = tximport(
            [salmon_file],
            "salmon",
            data_directory / "fabry_disease" / "transcript_gene_mapping_human.tsv",
            return_transcript_data=True,
            ignore_transcript_version=True,
            ignore_after_bar=True,
            output_type="xarray",
            counts_from_abundance=counts_from_abundance,  # type: ignore
        )

        assert isinstance(result, xr.Dataset), "The result is not an xarray Dataset."

        df_result = pd.DataFrame(
            result["counts"].data,
            index=result.coords["transcript_id"],
            columns=result.coords["file_path"],
        ).sort_index()

        # Load in the comparison data generated with the R package tximport
        if counts_from_abundance == "dtu_scaled_tpm":
            df_result_tximport = pd.read_csv(
                data_directory / "salmon" / "counts_tximport_dtuScaledTPM.csv",
                index_col=0,
                header=0,
            ).sort_index()
        elif counts_from_abundance == "scaled_tpm":
            df_result_tximport = pd.read_csv(
                data_directory / "salmon" / "counts_tximport_scaledTPM.csv",
                index_col=0,
                header=0,
            ).sort_index()

        # Remove the transcript version from the index since tximport does not remove it
        df_result_tximport.index = df_result_tximport.index.str.split(".").str[0]

        # Since tximport does not use the file path as column, use the columns from the comparison data
        df_result.columns = df_result_tximport.columns

        # Check that the data is the same
        pd.testing.assert_frame_equal(df_result, df_result_tximport)


def test_correctness_gene_level(
    rsem_files: Path,
) -> None:
    """Test that gene level input data is handled correctly.

    Args:
        rsem_files (Path): Path to the RSEM quantification files.
    """
    data_directory = rsem_files[0].parent

    result = tximport(
        rsem_files,
        "rsem",
        data_directory / "fabry_disease" / "transcript_gene_mapping_human.tsv",
        gene_level=True,
        ignore_transcript_version=True,
        ignore_after_bar=True,
        output_type="xarray",
    )

    assert isinstance(result, xr.Dataset), "The result is not an xarray Dataset."

    df_result = pd.DataFrame(
        result["counts"].data,
        index=result.coords["gene_id"],
        columns=result.coords["file_path"],
    ).sort_index()

    # Load in the comparison data generated with the R package tximport
    df_result_tximport = pd.read_csv(
        data_directory / "rsem" / "counts_tximport_no.csv",
        index_col=0,
        header=0,
    ).sort_index()

    # Remove the transcript version from the index since tximport does not remove it
    df_result_tximport.index = df_result_tximport.index.str.split(".").str[0]

    # Since tximport does not use the file path as column, use the columns from the comparison data
    df_result.columns = df_result_tximport.columns

    # Check that the data is the same
    pd.testing.assert_frame_equal(df_result, df_result_tximport)


def test_correctness_inferential_replicates(
    fabry_disease_files: List[Path],
) -> None:
    """Test importing salmon and kallisto quantification files with inferential replicates.

    Args:
        fabry_disease_files (List[Path]): The paths to the quantification files.
    """
    fabry_directory = fabry_disease_files[0].parent

    for data_type in ["kallisto", "salmon"]:
        for return_transcript_data in [True, False]:
            result = tximport(
                fabry_disease_files,
                data_type,  # type: ignore
                fabry_directory / "transcript_gene_mapping_human.tsv",
                return_transcript_data=return_transcript_data,
                inferential_replicates=True,
                inferential_replicate_transformer=lambda x: np.median(x, axis=1),
                inferential_replicate_variance=(data_type == "kallisto"),
                ignore_transcript_version=True,
                ignore_after_bar=True,
                output_type="xarray",
                counts_from_abundance=None,  # type: ignore
            )

            assert isinstance(result, xr.Dataset), "The result is not an xarray Dataset."

            df_result = pd.DataFrame(
                result["counts"].data,
                index=(result.coords["gene_id"] if not return_transcript_data else result.coords["transcript_id"]),
                columns=result.coords["file_path"],
            ).sort_index()
            df_result.index = df_result.index.str.split(".").str[0]

            # Load in the comparison data generated with the R package tximport
            if return_transcript_data:
                prefix = "_transcripts"
            else:
                prefix = ""

            df_result_tximport = pd.read_csv(
                fabry_directory / f"counts_tximport_bootstrap{prefix}_{data_type}.csv",
                index_col=0,
                header=0,
            )
            df_result_tximport.index = df_result_tximport.index.str.split(".").str[0]
            df_result_tximport = df_result_tximport.sort_index()

            # Since tximport does not use the file path as column, use the columns from the comparison data
            df_result.columns = df_result_tximport.columns

            # Check that the data is the same
            pd.testing.assert_frame_equal(df_result, df_result_tximport)
