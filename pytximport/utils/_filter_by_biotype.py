from logging import log
from pathlib import Path
from typing import List, Optional, Union

import anndata as ad
import numpy as np
import pandas as pd
import xarray as xr


def filter_by_biotype(
    transcript_data: Union[xr.Dataset, ad.AnnData],
    transcript_gene_map: Optional[Union[pd.DataFrame, Path, str]] = None,
    biotype_filter: Optional[List[str]] = None,
    id_column: str = "transcript_id",
    recalculate_abundance: bool = False,
) -> Union[xr.Dataset, ad.AnnData]:
    """Filter the transcripts by biotype.

    This function filters the transcripts by biotype. The biotype is assumed to be present in the transcript_id
    separated by a bar. The biotype is checked against the biotype_filter and the transcripts that match the biotype
    are kept. This function is provided mainly for internal use if `biotype_filter` is provided to the main function.

    Args:
        transcript_data (Union[xr.Dataset, ad.AnnData]): The expression data.
        transcript_gene_map (Union[pd.DataFrame, Path, str], optional): The mapping from transcript to gene with the
            `gene_biotype` column. If None, the biotype is assumed to be present in the id_column.
            Defaults to None.
        biotype_filter (List[str]): The biotypes to keep. Defaults to None.
        id_column (str, optional): The column name for the transcript/gene ID. Defaults to "transcript_id".
        recalculate_abundance (bool, optional): Whether to recalculate the abundance after filtering. This converts the
            abundance to TPM of the remaining transcripts but has implications for how the abundance can be used
            statistically. Defaults to False.

    Returns:
        Union[xr.Dataset, ad.AnnData]: The expression data filtered by biotype.
    """
    if isinstance(transcript_data, (str)):
        transcript_data = Path(transcript_data)

    if isinstance(transcript_data, ad.AnnData):
        transcript_ids = transcript_data.var_names.to_list()
    elif isinstance(transcript_data, xr.Dataset):
        transcript_ids = transcript_data.coords[id_column].values
    else:
        raise ValueError("The data type is not supported.")

    if biotype_filter is None:
        log(25, "No biotype filter provided. Skipping the biotype filtering.")
        return transcript_data

    if isinstance(transcript_gene_map, (str, Path)):
        transcript_gene_map = pd.read_csv(transcript_gene_map, header=0, index_col=None)

    if transcript_gene_map is None:
        # We can filter the transcripts by biotype if the biotype is present in the transcript_id since some
        # quantification tools include the biotype in the transcript_id
        # This only works if the transcript_id contains the biotype as one of the bar-separated fields and is not the
        # best way to filter the transcripts by biotype
        assert any("|" in transcript_id for transcript_id in transcript_ids), (
            "The transcript_id column does not contain the biotype. Please use the `pytximport.utils.filter_by_biotype`"
            " function with the `transcript_gene_map` argument instead. This function can be called after the"
            " `tximport` function using the resulting AnnData object or xarray Dataset."
        )

        transcript_id_fields = [transcript_id.split("|") for transcript_id in transcript_ids]

        transcript_keep_boolean = np.zeros_like(transcript_ids, dtype=bool)

        for biotype in biotype_filter:
            transcript_keep_boolean = np.logical_or(
                transcript_keep_boolean,
                [(biotype in transcript_id_field) for transcript_id_field in transcript_id_fields],
            )
    else:
        assert id_column in transcript_gene_map.columns, f"The {id_column} column is not present in the mapping file."
        assert (
            "gene_biotype" in transcript_gene_map.columns
        ), "The gene_biotype column is not present in the mapping file."

        transcript_id_fields = [transcript_id.split("|") for transcript_id in transcript_gene_map[id_column].values]

        # Filter the transcript-to-gene mapping by the biotype
        transcripts_to_keep = transcript_gene_map[transcript_gene_map["gene_biotype"].isin(biotype_filter)][
            [id_column]
        ].values

        # Filter the transcript data by the transcripts to keep
        transcript_keep_boolean = np.isin(transcript_ids, transcripts_to_keep)

    # Check that at least one transcript is protein-coding
    assert any(transcript_keep_boolean), "No transcript/gene with the desired biotype are present in the data."

    if isinstance(transcript_data, xr.Dataset):
        # Calculate the total abundance before filtering
        total_abundance = transcript_data["abundance"].sum(axis=0)
        transcript_data = transcript_data.isel(
            indexers={id_column: transcript_keep_boolean},
            drop=True,
        )

    elif isinstance(transcript_data, ad.AnnData):
        # Calculate the total abundance before filtering
        total_abundance = transcript_data.obsm["abundance"].sum(axis=1)
        transcript_data = transcript_data[:, transcript_keep_boolean]

    log(
        25,
        f"Removed {len(transcript_keep_boolean) - sum(transcript_keep_boolean)} transcripts with other biotypes.",
    )

    # Recalculate the abundance for each sample
    if recalculate_abundance:
        log(25, "Recalculating the abundance after filtering by biotype.")

        if isinstance(transcript_data, ad.AnnData):
            new_abundance = transcript_data.obsm["abundance"].sum(axis=1)
            ratio = total_abundance / new_abundance
            transcript_data.obsm["abundance"] = (transcript_data.obsm["abundance"].T * ratio).T
        elif isinstance(transcript_data, xr.Dataset):
            new_abundance = transcript_data["abundance"].sum(axis=0)
            ratio = total_abundance / new_abundance
            transcript_data["abundance"] = (transcript_data["abundance"].T * ratio).T

    return transcript_data
