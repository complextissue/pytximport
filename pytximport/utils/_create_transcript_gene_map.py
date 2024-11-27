import re
from logging import log, warning
from pathlib import Path
from typing import Any, Dict, List, Literal, Union

import numpy as np
import pandas as pd


def create_transcript_gene_map(
    species: Literal["human", "mouse"] = "human",
    host: str = "http://www.ensembl.org",
    source_field: Literal["ensembl_transcript_id", "external_transcript_name"] = "ensembl_transcript_id",
    target_field: Union[
        Literal["ensembl_gene_id", "external_gene_name", "external_transcript_name", "gene_biotype"],
        List[Literal["ensembl_gene_id", "external_gene_name", "external_transcript_name", "gene_biotype"]],
    ] = "ensembl_gene_id",
    rename_columns: bool = True,
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

        # or get multiple fields
        transcript_gene_map = create_transcript_gene_map(
            species="mouse",
            target_field=["external_gene_name", "gene_biotype"],
        )

    Args:
        species (Literal["human", "mouse"], optional): The species to use. Defaults to "human".
        host (str, optional): The host to use. Defaults to "http://www.ensembl.org".
        source_field (Literal["ensembl_transcript_id", "external_transcript_name"], optional): The identifier to get for
            each transcript id. Defaults to "ensembl_transcript_id".
        target_field (Union[Literal["ensembl_gene_id", "external_gene_name", "external_transcript_name", "gene_biotype"]
            , List[Literal["ensembl_gene_id", "external_gene_name", "external_transcript_name", "gene_biotype"]]],
            optional): The corresponding identifier to get for each transcript. Defaults to "ensembl_gene_id".
        rename_columns (bool, optional): Whether to rename `ensembl_transcript_id` to `transcript_id`, `ensembl_gene_id`
            to `gene_id`, `external_gene_name` to `gene_name` if the gene id is also present or `gene_id` if no other
            gene id is present, and `external_transcript_name` to `transcript_name`. Defaults to True.

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

    columns: List[str] = [source_field]

    if isinstance(target_field, list):
        columns = [*columns, *target_field]
    else:
        columns = [*columns, target_field]

    transcript_gene_map: pd.DataFrame = dataset.query(attributes=columns)
    transcript_gene_map.columns = pd.Index(columns)

    if rename_columns:
        transcript_gene_map.rename(
            columns={
                "ensembl_transcript_id": "transcript_id",
                "external_transcript_name": "transcript_name",
                "ensembl_gene_id": "gene_id",
                "external_gene_name": (
                    "gene_name"
                    if target_field == "ensembl_gene_id"
                    or (isinstance(target_field, list) and "ensembl_gene_id" in target_field)
                    else "gene_id"
                ),
            },
            inplace=True,
        )

    transcript_gene_map.dropna(inplace=True)
    transcript_gene_map.drop_duplicates(inplace=True)
    transcript_gene_map.reset_index(drop=True, inplace=True)

    return transcript_gene_map


def create_transcript_gene_map_from_annotation(
    file_path: Union[str, Path],
    source_field: Literal["transcript_id", "transcript_name"] = "transcript_id",
    target_field: Union[
        Literal["gene_id", "gene_name", "gene_biotype", "transcript_name"],
        List[Literal["gene_id", "gene_name", "gene_biotype", "transcript_name"]],
    ] = "gene_id",
    use_transcript_name_as_replacement_id: bool = True,
    use_gene_name_as_replacement_id: bool = True,
    chunk_size: int = 100000,
    **kwargs: Dict[str, Any],
) -> pd.DataFrame:
    """Create a mapping from transcript ids to gene ids using a GTF annotation file.

    Basic example:

    .. code-block:: python

            from pytximport.utils import create_transcript_gene_map_from_annotation

            # Create a mapping from transcript ids to gene names
            transcript_gene_map = create_transcript_gene_map_from_annotation(
                "path/to/annotation.gtf",
                target_field="gene_name",
            )

            # Create a mapping from transcript ids to transcript names and include the gene biotype
            transcript_gene_map = create_transcript_gene_map_from_annotation(
                "path/to/annotation.gtf",
                target_field=["transcript_name", "gene_biotype"],
            )

    Args:
        file_path (Union[str, Path]): The path to the GTF annotation file.
        source_field (Literal["transcript_id", "transcript_name"], optional): The identifier to get for each transcript
            id. Defaults to "transcript_id".
        target_field (Union[ Literal["gene_id", "gene_name", "gene_biotype"], List[Literal["gene_id", "gene_name",
            "gene_biotype"]], optional): The corresponding identifier(s) to get for each transcript.
            Defaults to "gene_id".
        chunk_size (int, optional): The number of lines to read at a time. Defaults to 100000.
        use_transcript_name_as_replacement_id (bool, optional): Whether to use the transcript name as the transcript id
            if the transcript id is missing. Defaults to True.
        use_gene_name_as_replacement_id (bool, optional): Whether to use the gene name as the gene id if the gene id is
            missing. Defaults to True.
        keep_biotype (bool, optional): Whether to keep the gene_biotype column. Defaults to False.

    Returns:
        pd.DataFrame: The mapping from transcript ids to gene ids.
    """
    assert source_field != target_field, "The source_field and target_field must be different."

    transcript_gene_map = pd.DataFrame(columns=["transcript_id", "gene_id", "gene_name", "gene_biotype"])

    if "field" in kwargs:
        warning("The field argument is deprecated. Please use the source_field and target_field arguments instead.")

    if "keep_biotype" in kwargs and kwargs["keep_biotype"]:
        warning("The keep_biotype argument is deprecated. Please use the target_field argument with a list instead.")
        if target_field != "gene_biotype" and not (isinstance(target_field, list) and "gene_biotype" in target_field):
            target_field = (
                [*target_field, "gene_biotype"] if isinstance(target_field, list) else [target_field, "gene_biotype"]
            )

    if not isinstance(file_path, Path):
        file_path = Path(file_path)

    if not Path(file_path).exists():
        raise FileNotFoundError(f"The file {file_path} does not exist.")

    for chunk in pd.read_csv(
        file_path,
        sep="\t",
        chunksize=chunk_size,
        header=None,
        comment="#",
        compression=("gzip" if file_path.suffix == ".gz" else None),
        engine="c",
    ):
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

    # If only the transcript_name is present, we can drop the id and rename the transcript_name to transcript_id
    if source_field == "transcript_name" and use_transcript_name_as_replacement_id:
        transcript_gene_map.drop("transcript_id", axis=1, inplace=True)
        transcript_gene_map.rename(columns={"transcript_name": "transcript_id"}, inplace=True)

        source_field = "transcript_id"

    # If only the gene_name is present, we can drop the gene_id and rename the gene_name to gene_id
    if (
        (target_field == "gene_name" or (isinstance(target_field, list) and "gene_name" in target_field))
        and not (target_field == "gene_id" or (isinstance(target_field, list) and "gene_id" in target_field))
        and use_gene_name_as_replacement_id
    ):
        log(
            25,
            (
                "No gene_id target field was provided. Renaming gene_name to gene_id. "
                "You can disable this behavior by setting use_gene_name_as_replacement_id to False."
            ),
        )

        transcript_gene_map.drop("gene_id", axis=1, inplace=True)
        transcript_gene_map.rename(columns={"gene_name": "gene_id"}, inplace=True)

        if isinstance(target_field, list):
            target_field = [field if field != "gene_name" else "gene_id" for field in target_field]
        else:
            target_field = "gene_id" if target_field == "gene_name" else target_field

    fields_to_keep = [source_field, *target_field] if isinstance(target_field, list) else [source_field, target_field]

    transcript_gene_map = transcript_gene_map[fields_to_keep]
    transcript_gene_map.replace("", np.nan, inplace=True)
    transcript_gene_map.dropna(inplace=True)

    if source_field == "transcript_id" and (
        target_field == "gene_id" or (isinstance(target_field, list) and "gene_id" in target_field)
    ):
        transcript_gene_map.drop_duplicates(subset=["transcript_id", "gene_id"], inplace=True)

    transcript_gene_map.reset_index(drop=True, inplace=True)

    return transcript_gene_map
