"""Test importing salmon quantification files."""

from pathlib import Path
from typing import List

import anndata as ad
import pandas as pd
import xarray as xr

from pytximport import tximport
from pytximport.utils import replace_transcript_ids_with_names


def test_salmon_transcript_level(
    salmon_file: Path,
    transcript_name_mapping_human: pd.DataFrame,
) -> None:
    """Test importing a salmon quantification file.

    Args:
        salmon_file (Path): Path to the salmon quantification file.
    """
    result = tximport(
        [salmon_file],
        "salmon",
        return_transcript_data=True,
        output_type="anndata",
        ignore_transcript_version=True,
        ignore_after_bar=True,
        counts_from_abundance="scaled_tpm",
    )

    assert isinstance(result, ad.AnnData)
    counts = result.X

    # check that the counts.data are all positive
    assert (counts >= 0).all()

    # replace the transcript ids with the transcript names
    result = replace_transcript_ids_with_names(result, transcript_name_mapping_human)

    # check that the result is an AnnData object
    assert isinstance(result, ad.AnnData)

    # check that the var names don't start with ENST
    assert result.var_names.str.startswith("ENST").sum() == 0

    # check that the var names are not nan or the string "nan" or empty
    assert (
        result.var_names.isna().sum() == 0
        and (result.var_names == "nan").sum() == 0
        and (result.var_names == "").sum() == 0
    )
