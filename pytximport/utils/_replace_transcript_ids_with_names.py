from pathlib import Path
from typing import Optional, Union

import anndata as ad
import pandas as pd
import xarray as xr

from ._remove_transcript_version import remove_transcript_version


def replace_transcript_ids_with_names(
    transcript_data: Union[ad.AnnData, xr.Dataset],
    transcript_name_map: Union[pd.DataFrame, Union[str, Path]],
) -> Union[ad.AnnData, xr.Dataset]:
    """Replace transcript IDs with transcript names.

    Args:
        transcript_data (Union[ad.AnnData, xr.Dataset]): The transcript-level expression data.
        transcript_name_map (Union[pd.DataFrame, Union[str, Path]]): The mapping from transcripts to
            names. Contains two columns: `transcript_id` and `transcript_name`.

    Returns:
        Union[ad.AnnData, xr.Dataset]: The transcript-level expression data with the transcript names.
    """
    # Read the transcript to gene mapping
    if isinstance(transcript_name_map, str) or isinstance(transcript_name_map, Path):
        transcript_name_map = pd.read_table(
            transcript_name_map,
            header=0,
            engine="c",
            usecols=["transcript_id", "transcript_name"],
            dtype=str,
        )
        transcript_name_map = transcript_name_map.drop_duplicates()

    # Check that transcript_id and transcript_name are present in the mapping
    assert "transcript_id" in transcript_name_map.columns, "The mapping does not contain a `transcript_id` column."
    assert "transcript_name" in transcript_name_map.columns, "The mapping does not contain a `transcript_name` column."

    # Check whether the transcript_data is an AnnData object and convert it to an xr.Dataset
    return_as_anndata = False
    if isinstance(transcript_data, ad.AnnData):
        return_as_anndata = True
        transcript_data = xr.Dataset(
            data_vars={
                "counts": xr.DataArray(transcript_data.X.T, dims=["transcript_id", "file"]),
                "length": xr.DataArray(transcript_data.obsm["length"].T, dims=["transcript_id", "file"]),
                "abundance": xr.DataArray(transcript_data.obsm["abundance"].T, dims=["transcript_id", "file"]),
            },
            coords={
                "transcript_id": transcript_data.var.index.values,
                "file_path": transcript_data.obs.index.values,
            },
        )

    # Remove the transcript version
    transcript_data, transcript_name_map, _ = remove_transcript_version(  # type: ignore
        transcript_data,
        transcript_name_map,
    )

    transcript_name_dict = transcript_name_map.set_index("transcript_id")["transcript_name"].to_dict()
    transcript_names = transcript_data["transcript_id"].to_series().map(transcript_name_dict).values

    # Replace the transcript_id with the transcript_name
    transcript_data.coords["transcript_id"] = transcript_names

    if return_as_anndata:
        return ad.AnnData(
            X=transcript_data["counts"].values.T,
            obs=pd.DataFrame(index=transcript_data.coords["file_path"].values),
            var=pd.DataFrame(index=transcript_data.coords["transcript_id"].values),
            obsm={
                "length": transcript_data["length"].values.T,
                "abundance": transcript_data["abundance"].values.T,
            },
        )

    return transcript_data
