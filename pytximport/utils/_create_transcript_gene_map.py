import re
from logging import warning
from pathlib import Path
from typing import Any, Dict, Literal, Union

import numpy as np
import pandas as pd


def create_transcript_gene_map(
    species: Literal["human", "mouse"] = "human",
    host: str = "http://www.ensembl.org",
    source_field: Literal["ensembl_transcript_id", "external_transcript_name"] = "ensembl_transcript_id",
    target_field: Literal["ensembl_gene_id", "external_gene_name", "external_transcript_name"] = "ensembl_gene_id",
    **kwargs: Dict[str, Any],
) -> pd.DataFrame:
    """Create a mapping from transcript ids to gene ids using the Ensembl Biomart.

    .. warning ::
        Choosing any `target_field` value other than `ensembl_gene_id` may not result in a full transcript to gene map
        since not all transcripts may have the respective variable. While this does not typically affect well defined
        transcripts, be aware of this possible source of bias.

    Basic example:

    .. code-block:: python

        from pytximport.utils import create_transcript_gene_map

        transcript_gene_map = create_transcript_gene_map(
            species="human",
            host="https://may2024.archive.ensembl.org/", # Use a specific Ensembl release
            target_field="external_gene_name",
        )

    Args:
        species (Literal["human", "mouse"], optional): The species to use. Defaults to "human".
        host (str, optional): The host to use. Defaults to "http://www.ensembl.org".
        source_field (Literal["ensembl_transcript_id", "external_transcript_name"], optional): The identifier to get for
            each transcript id. Defaults to "ensembl_transcript_id".
        target_field (Literal["ensembl_gene_id", "external_gene_name", "external_transcript_name"], optional): The
            corresponding identifier to get for each transcript. Defaults to "ensembl_gene_id".

    Returns:
        pd.DataFrame: The mapping from transcript ids to gene ids.
    """
    from pybiomart import Dataset

    if "field" in kwargs:
        warning("The field argument is deprecated. Please use the source_field and target_field arguments instead.")

    if species == "human":
        dataset = Dataset(name="hsapiens_gene_ensembl", host=host)
    elif species == "mouse":
        dataset = Dataset(name="mmusculus_gene_ensembl", host=host)

    transcript_gene_map: pd.DataFrame = dataset.query(attributes=[source_field, target_field])
    transcript_gene_map.columns = pd.Index(
        [
            "transcript_id",
            ("gene_id" if target_field != "external_transcript_name" else "transcript_name"),
        ]
    )

    transcript_gene_map.dropna(inplace=True)
    transcript_gene_map.drop_duplicates(inplace=True)
    transcript_gene_map.reset_index(drop=True, inplace=True)

    return transcript_gene_map


def create_transcript_gene_map_from_annotation(
    file_path: Union[str, Path],
    source_field: Literal["transcript_id", "transcript_name"] = "transcript_id",
    target_field: Literal["gene_id", "gene_name"] = "gene_id",
    chunk_size: int = 100000,
    keep_biotype: bool = False,
    **kwargs: Dict[str, Any],
) -> pd.DataFrame:
    """Create a mapping from transcript ids to gene ids using a GTF annotation file.

    Args:
        file_path (Union[str, Path]): The path to the GTF annotation file.
        field (Literal["gene_id", "gene_name"], optional): The identifier to get for each transcript id.
            Defaults to "gene_id".
        chunk_size (int, optional): The number of lines to read at a time. Defaults to 100000.
        keep_biotype (bool, optional): Whether to keep the gene_biotype column. Defaults to False.

    Returns:
        pd.DataFrame: The mapping from transcript ids to gene ids.
    """
    transcript_gene_map = pd.DataFrame(columns=["transcript_id", "gene_id", "gene_name", "gene_biotype"])

    if "field" in kwargs:
        warning("The field argument is deprecated. Please use the source_field and target_field arguments instead.")

    for chunk in pd.read_csv(file_path, sep="\t", chunksize=chunk_size, header=None, comment="#"):
        # See: https://www.ensembl.org/info/website/upload/gff.html
        chunk.columns = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]

        # Each attribute line looks like this:
        # gene_id ""; transcript_id ""; gene_name ""; gene_source ""; gene_biotype "";
        # transcript_name ""; transcript_source "";
        # We are only interested in the gene_id, gene_name, transcript_id, transcript_name and gene_biotype
        attribute_columns = [
            "transcript_id",
            "transcript_name",
            "gene_id",
            "gene_name",
            "gene_biotype",
        ]
        for column in attribute_columns:
            chunk[column] = chunk["attribute"].apply(
                # Adopted from https://gist.github.com/rf-santos/22f521c62ca2f85ac9582bf0d91e4054
                lambda x: (re.findall(rf'{column} "([^"]*)"', x)[0] if rf'{column} "' in x else "")
            )

        chunk.drop("attribute", axis=1, inplace=True)

        transcript_gene_map = pd.concat(
            [
                transcript_gene_map,
                chunk[["transcript_id", "transcript_name", "gene_id", "gene_name", "gene_biotype"]],
            ]
        )

    # Replace the gene_name with the gene_id where the gene_name is ""
    transcript_gene_map["gene_name"] = np.where(
        transcript_gene_map["gene_name"] == "",
        transcript_gene_map["gene_id"],
        transcript_gene_map["gene_name"],
    )

    if source_field == "transcript_name":
        transcript_gene_map.drop("transcript_id", axis=1, inplace=True)
        transcript_gene_map.rename(columns={"transcript_name": "transcript_id"}, inplace=True)

    if target_field == "gene_name":
        transcript_gene_map.drop("gene_id", axis=1, inplace=True)
        transcript_gene_map.rename(columns={"gene_name": "gene_id"}, inplace=True)

    fields_to_keep = ["transcript_id", "gene_id"]

    if keep_biotype and "gene_biotype" in transcript_gene_map.columns:
        fields_to_keep.append("gene_biotype")

    transcript_gene_map = transcript_gene_map[fields_to_keep]

    transcript_gene_map[["gene_id", "transcript_id"]] = transcript_gene_map[["gene_id", "transcript_id"]].replace(
        "", np.nan
    )
    transcript_gene_map.dropna(inplace=True)
    transcript_gene_map.drop_duplicates(subset=["gene_id", "transcript_id"], inplace=True)
    transcript_gene_map.reset_index(drop=True, inplace=True)

    return transcript_gene_map
