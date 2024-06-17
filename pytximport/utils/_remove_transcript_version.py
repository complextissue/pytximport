from typing import List, Optional, Tuple

import pandas as pd
import xarray as xr


def remove_transcript_version(
    transcript_data: xr.Dataset,
    transcript_target_map: Optional[pd.DataFrame] = None,
    transcript_ids: Optional[List[str]] = None,
) -> Tuple[xr.Dataset, pd.DataFrame, List[str]]:
    """Remove the transcript version from the transcript data and the transcript target map.

    Args:
        transcript_data (xr.Dataset): The transcript data.
        transcript_target_map (Optional[pd.DataFrame], optional): The transcript target map. Defaults to None.
        transcript_ids (Optional[List[str]], optional): The transcript ids. Defaults to None.

    Returns:
        Tuple[xr.Dataset, pd.DataFrame, List[str]]: The transcript data, the transcript target map, and the transcript
        ids.
    """
    # ignore the transcript version in both the data and the transcript gene map
    if transcript_ids is None:
        transcript_ids = transcript_data.coords["transcript_id"].values  # type: ignore

    if transcript_target_map is not None:
        transcript_target_map["transcript_id"] = transcript_target_map["transcript_id"].str.split(".").str[0]

    transcript_ids = [transcript_id.split(".")[0] for transcript_id in transcript_ids]  # type: ignore
    transcript_data.coords["transcript_id"] = transcript_ids

    return transcript_data, transcript_target_map, transcript_ids
