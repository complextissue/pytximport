from logging import warning
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple, Union

import numpy as np
import xarray as xr
from tqdm import tqdm


def summarize_rsem_gene_data(
    file_paths: Union[List[str], List[Path]],
    importer: Callable,
    importer_kwargs: Dict[str, Any],
    existence_optional: bool = False,
) -> Tuple[Optional[xr.Dataset], List[int]]:
    """Summarize gene-level RSEM quantification files.

    Args:
        file_paths (Union[List[str], List[Path]]): The paths to the quantification files.
        importer (Callable): The importer function to read the quantification files.
        importer_kwargs (Dict[str, Any]): The keyword arguments for the importer function.
        existence_optional (bool, optional): Whether the files are optional. Defaults to False.

    Return:
        xr.Dataset: The gene-level expression data.
    """
    transcript_data: Optional[xr.Dataset] = None
    file_paths_missing_idx: List[int] = []

    for file_idx, file_path in tqdm(enumerate(file_paths), desc="Reading quantification files"):
        try:
            transcript_data_sample = importer(file_path, **importer_kwargs)
        except ImportError as exception:
            if existence_optional:
                file_paths_missing_idx.append(file_idx)
                warning(f"Could not import the file: {file_path}")
                continue
            else:
                raise exception

        if file_idx == 0 or transcript_data is None:
            # Since the files are parsed with the generic importer, the gene id column is called "transcript_id"
            empty_array = np.zeros((len(transcript_data_sample["transcript_id"]), len(file_paths)))

            abundance = xr.DataArray(data=empty_array.copy(), dims=["gene_id", "file"])
            counts = xr.DataArray(data=empty_array.copy(), dims=["gene_id", "file"])
            length = xr.DataArray(data=empty_array.copy(), dims=["gene_id", "file"])

            data_vars = {
                "abundance": abundance,
                "counts": counts,
                "length": length,
            }

            transcript_data = xr.Dataset(
                data_vars,
                coords={
                    "gene_id": transcript_data_sample["transcript_id"],
                    "file_path": list(file_paths),
                },
            )

        transcript_data["abundance"][:, file_idx] = transcript_data_sample["abundance"]
        transcript_data["counts"][:, file_idx] = transcript_data_sample["counts"]
        transcript_data["length"][:, file_idx] = transcript_data_sample["length"]

    return transcript_data, file_paths_missing_idx
