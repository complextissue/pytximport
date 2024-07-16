import re
from pathlib import Path
from typing import Literal, Union

import numpy as np
import pandas as pd


def create_transcript_to_gene_map(
    species: Literal["human", "mouse"] = "human",
    host: str = "http://www.ensembl.org",
    field: Literal["ensembl_gene_id", "external_gene_name", "external_transcript_name"] = "ensembl_gene_id",
) -> pd.DataFrame:
    """Create a mapping from transcript ids to gene ids using the Ensembl Biomart.

    Args:
        species (Literal["human", "mouse"], optional): The species to use. Defaults to "human".
        host (str, optional): The host to use. Defaults to "http://www.ensembl.org".
        field (Literal["ensembl_gene_id", "external_gene_name", "external_transcript_name"], optional): The
            identifier to get for each transcript id. Defaults to "external_gene_name".

    Returns:
        pd.DataFrame: The mapping from transcript ids to gene ids.
    """
    from pybiomart import Dataset

    if species == "human":
        dataset = Dataset(name="hsapiens_gene_ensembl", host=host)
    elif species == "mouse":
        dataset = Dataset(name="mmusculus_gene_ensembl", host=host)

    transcript_gene_map = dataset.query(attributes=["ensembl_transcript_id", field])
    transcript_gene_map.columns = [
        "transcript_id",
        ("gene_id" if field != "external_transcript_name" else "transcript_name"),
    ]

    transcript_gene_map.dropna(inplace=True)
    transcript_gene_map.drop_duplicates(inplace=True)
    transcript_gene_map.reset_index(drop=True, inplace=True)

    return transcript_gene_map


def create_transcript_to_gene_map_from_gtf_annotation(
    file_path: Union[str, Path],
    field: Literal["gene_id", "gene_name"] = "gene_id",
    chunk_size: int = 100000,
    keep_gene_name: bool = True,
    keep_biotype: bool = False,
) -> pd.DataFrame:
    """Create a mapping from transcript ids to gene ids using a GTF annotation file.

    Args:
        file_path (Union[str, Path]): The path to the GTF annotation file.
        field (Literal["gene_id", "gene_name"], optional): The identifier to get for each transcript id.
            Defaults to "gene_id".
        chunk_size (int, optional): The number of lines to read at a time. Defaults to 100000.
        keep_gene_name (bool, optional): Whether to keep the gene_name column when field is "gene_id".
            Defaults to True.
        keep_biotype (bool, optional): Whether to keep the gene_biotype column. Defaults to False.

    Returns:
        pd.DataFrame: The mapping from transcript ids to gene ids.
    """
    transcript_gene_map = pd.DataFrame(columns=["transcript_id", "gene_id", "gene_name", "gene_biotype"])

    for chunk in pd.read_csv(file_path, sep="\t", chunksize=chunk_size, header=None, comment="#"):
        # see: https://www.ensembl.org/info/website/upload/gff.html
        chunk.columns = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]

        # each attribute line looks like this:
        # gene_id ""; transcript_id ""; gene_name ""; gene_source ""; gene_biotype "";
        # transcript_name ""; transcript_source "";
        # however, we are only interested in the gene_id, gene_name, and transcript_id
        attribute_columns = [
            "transcript_id",
            "gene_id",
            "gene_name",
            "gene_biotype",
        ]
        for column in attribute_columns:
            chunk[column] = chunk["attribute"].apply(
                # adopted from https://gist.github.com/rf-santos/22f521c62ca2f85ac9582bf0d91e4054
                lambda x: (re.findall(rf'{column} "([^"]*)"', x)[0] if rf'{column} "' in x else "")
            )

        chunk.drop("attribute", axis=1, inplace=True)

        transcript_gene_map = pd.concat(
            [
                transcript_gene_map,
                chunk[["transcript_id", "gene_id", "gene_name", "gene_biotype"]],
            ]
        )

    # replace the gene_name with the gene_id where the gene_name is ""
    transcript_gene_map["gene_name"] = np.where(
        transcript_gene_map["gene_name"] == "",
        transcript_gene_map["gene_id"],
        transcript_gene_map["gene_name"],
    )

    if field == "gene_name":
        transcript_gene_map.drop("gene_id", axis=1, inplace=True)
        transcript_gene_map.rename(columns={"gene_name": "gene_id"}, inplace=True)

    if not keep_gene_name and "gene_name" in transcript_gene_map.columns:
        transcript_gene_map.drop("gene_name", axis=1, inplace=True)

    if not keep_biotype and "gene_biotype" in transcript_gene_map.columns:
        transcript_gene_map.drop("gene_biotype", axis=1, inplace=True)

    transcript_gene_map[["gene_id", "transcript_id"]] = transcript_gene_map[["gene_id", "transcript_id"]].replace(
        "", np.nan
    )
    transcript_gene_map.dropna(inplace=True)
    transcript_gene_map.drop_duplicates(subset=["gene_id", "transcript_id"], inplace=True)
    transcript_gene_map.reset_index(drop=True, inplace=True)

    return transcript_gene_map
