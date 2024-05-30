from logging import log, warning
from pathlib import Path
from time import time
from typing import Callable, List, Literal, Optional, Union

import anndata as ad
import numpy as np
import pandas as pd
import xarray as xr
from tqdm import tqdm

from ..importers import read_kallisto, read_salmon, read_tsv
from ..utils import convert_abundance_to_counts, convert_transcripts_to_genes


def tximport(
    file_paths: Union[List[str], List[Path]],
    data_type: Literal["kallisto", "salmon", "sailfish", "oarfish", "piscem", "stringtie", "tsv"],
    transcript_gene_map: Optional[Union[pd.DataFrame, Union[str, Path]]] = None,
    # the order of transcript_gene_map and counts_from_abundance is swapped, since the former is required for most cases
    counts_from_abundance: Optional[Literal["scaled_tpm", "length_scaled_tpm"]] = None,
    # no equivalent to txIn, as RSEM is not currently supported and all other data are transcripts
    return_transcript_data: bool = False,
    # inferential replicates are not supported but the argument is kept for compatibility
    inferential_replicates: bool = False,
    inferential_replicate_transformer: Optional[Callable] = None,
    inferential_replicate_variance_reduce: bool = False,
    ignore_transcript_version: bool = True,
    ignore_after_bar: bool = True,
    id_column: Optional[str] = None,
    counts_column: Optional[str] = None,
    length_column: Optional[str] = None,
    abundance_column: Optional[str] = None,
    custom_importer: Optional[Callable] = None,
    existence_optional: bool = False,
    sparse: bool = False,
    sparse_threshold: Optional[float] = None,
    read_length: Optional[int] = None,
    # arguments exclusive to the pytximport implementation
    output_type: Literal["xarray", "anndata"] = "xarray",
    output_format: Literal["csv", "h5ad"] = "csv",
    save_path: Optional[Union[str, Path]] = None,
    return_data: bool = True,
    biotype_filter: Optional[List[str]] = None,
) -> Union[xr.Dataset, ad.AnnData, None]:
    """Import transcript-level quantification files and convert them to gene-level expression.

    An alternative to the `tximport` R package. `pytximport` is a Python package that imports transcript-level
    quantification files and converts them to gene-level expression. Currently, only `kallisto` quantification files are
    supported.

    .. note::
        Please also cite the original `tximport` R package when using `pytximport`.
        DOI: https://doi.org/10.12688/f1000research.7563.1

    .. warning::
        The `pytximport` package is still in development and currently does not make use of the inferential replicates,
        unlike the `tximport` R package. Default parameters, argument names and return types may also differ between the
        implementations.

    Args:
        file_paths (List[Union[str, Path]]): The paths to the quantification files.
        data_type (Literal["kallisto", "salmon"], optional): The type of quantification file.
        transcript_gene_map (pd.DataFrame): The mapping from transcripts to genes. Contains two columns: `transcript`
            and `gene`.
        counts_from_abundance (Optional[Literal["scaled_tpm", "length_scaled_tpm"]], optional): The type of counts to
            convert to. Defaults to "length_scaled_tpm".
        return_transcript_data (bool, optional): Whether to return the transcript-level expression. Defaults to False.
        inferential_replicates (bool, optional): Whether to drop the inferential replicates. Note, that this
            implementation currently has no support for inferential replicates. Defaults to False.
        inferential_replicate_transformer (Optional[Callable], optional): A custom function to transform the inferential
            replicates. Currently unsupported. Defaults to None.
        inferential_replicate_variance_reduce (bool, optional): Whether to replace the inferential replicates with the
            sample variance for each transcript. Currently unsupported. Defaults to False.
        ignore_transcript_version (bool, optional): Whether to ignore the transcript version. Defaults to True.
        ignore_after_bar (bool, optional): Whether to split the transcript id after the bar character (`|`).
            Defaults to True.
        id_column (Optional[str], optional): The column name for the transcript id. Defaults to None.
        counts_column (Optional[str], optional): The column name for the counts. Defaults to None.
        length_column (Optional[str], optional): The column name for the length. Defaults to None.
        abundance_column (Optional[str], optional): The column name for the abundance. Defaults to None.
        custom_importer (Optional[Callable], optional): A custom importer function. Defaults to None.
        existence_optional (bool, optional): Whether the existence of the files is optional. Defaults to False.
        sparse (bool, optional): Whether to use sparse matrices. Currenlty, sparse input is not supported.
            Defaults to False.
        sparse_threshold (Optional[float], optional): The threshold for the sparse matrix. Currently, sparse input is
            not supported. Defaults to None.
        read_length (Optional[int], optional): The read length for the stringtie quantification. Defaults to None.
        output_type (Literal["xarray", "anndata"], optional): The type of output. Defaults to "xarray".
        output_format (Literal["csv", "h5ad"], optional): The type of output file. Defaults to "csv".
        save_path (Optional[Union[str, Path]], optional): The path to save the gene-level expression. Defaults to None.
        return_data (bool, optional): Whether to return the gene-level expression. Defaults to True.
        biotype_filter (List[str], optional): Filter the transcripts by biotype, including only those provided.
            Defaults to None.

    Returns:
        Union[xr.Dataset, ad.AnnData, None]: The estimated gene-level expression data if `return_data` is True.
    """
    # start a timer
    log(25, "Starting the import.")
    start_time = time()

    # needed for faster groupby operations
    xr.set_options(use_flox=True)

    if data_type == "sailfish":
        data_type = "salmon"

    assert data_type in [
        "kallisto",
        "salmon",
        "oarfish",
        "piscem",
        "stringtie",
        "tsv",
    ], "Only kallisto, salmon/sailfish, oarfish, piscem, stringtie, and tsv quantification files are supported."
    assert not sparse, "Currently, sparse matrices are not supported."
    assert not inferential_replicates, "Currently, inferential replicates are not supported."

    # read the transcript to gene mapping
    if isinstance(transcript_gene_map, str) or isinstance(transcript_gene_map, Path):
        transcript_gene_map = pd.read_table(transcript_gene_map, header=0)
        transcript_gene_map = transcript_gene_map.drop_duplicates()

    # assert that transcript_id and gene_id are present in the mapping
    if transcript_gene_map is not None:
        assert "transcript_id" in transcript_gene_map.columns, "The mapping does not contain a `transcript_id` column."
        assert "gene_id" in transcript_gene_map.columns, "The mapping does not contain a `gene_id` column."

    # assert that return_transcript_data is True if transcript_gene_map is None
    if transcript_gene_map is None:
        assert return_transcript_data, "A transcript to gene mapping must be provided when summarizing to genes."

    # create the importer kwargs based on the provided column names
    importer_kwargs = {
        "id_column": id_column,
        "counts_column": counts_column,
        "length_column": length_column,
        "abundance_column": abundance_column,
    }

    # list is required to create a copy so that the original dictionary can be changed
    for key, value in list(importer_kwargs.items()):
        if value is None:
            del importer_kwargs[key]

    # match the data type with the required importer
    importer: Callable
    if data_type == "kallisto":
        importer = read_kallisto
    elif data_type == "salmon":
        importer = read_salmon
    elif data_type == "oarfish":
        # use the default column names if not provided
        id_column = "tname" if id_column is None else id_column
        counts_column = "num_reads" if counts_column is None else counts_column
        length_column = "len" if length_column is None else length_column
        abundance_column = "num_reads" if abundance_column is None else abundance_column

        importer = read_tsv
    elif data_type == "piscem":
        id_column = "Name" if id_column is None else id_column
        counts_column = "NumReads" if counts_column is None else counts_column
        length_column = "EffectiveLength" if length_column is None else length_column

        importer = read_tsv
    elif data_type == "stringtie":
        id_column = "t_name" if id_column is None else id_column
        counts_column = "cov" if counts_column is None else counts_column
        length_column = "length" if length_column is None else length_column
        abundance_column = "FPKM" if abundance_column is None else abundance_column

        if read_length is None:
            raise ValueError("The read_length must be provided for stringtie quantification files.")

        importer = read_tsv
    elif data_type == "tsv":
        # check that all of id_column counts_column length_column abundance_column are provided
        if any(
            [
                id_column is None,
                counts_column is None,
                length_column is None,
                abundance_column is None,
            ]
        ):
            raise ValueError(
                "The id_column, counts_column, length_column, and abundance_column must be provided for .tsv files."
            )

        importer = read_tsv
    else:
        raise ValueError("The input file type is not supported.")

    # overwrite the importer if a custom importer is provided
    if custom_importer is not None:
        importer = custom_importer

    transcript_data: Optional[xr.Dataset] = None
    file_paths_missing_idx: List[int] = []

    # iterate through the files
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

        # the check for transcript_data is needed in case the import of the first file fails
        if file_idx == 0 or transcript_data is None:
            empty_array = np.zeros((len(transcript_data_sample["transcript_id"]), len(file_paths)))

            abundance = xr.DataArray(data=empty_array.copy(), dims=["transcript_id", "file"])
            counts = xr.DataArray(data=empty_array.copy(), dims=["transcript_id", "file"])
            length = xr.DataArray(data=empty_array.copy(), dims=["transcript_id", "file"])

            transcript_data = xr.Dataset(
                {
                    "abundance": abundance,
                    "counts": counts,
                    "length": length,
                },
                coords={
                    "transcript_id": transcript_data_sample["transcript_id"],
                    "file_path": list(file_paths),
                },
            )

        # add the transcript-level expression to the array
        transcript_data["abundance"].loc[{"file": file_idx}] = transcript_data_sample["abundance"]
        transcript_data["counts"].loc[{"file": file_idx}] = transcript_data_sample["counts"]
        transcript_data["length"].loc[{"file": file_idx}] = transcript_data_sample["length"]

    if transcript_data is None:
        raise ImportError("No files could be imported.")

    if existence_optional and len(file_paths_missing_idx):
        if len(file_paths_missing_idx) == len(file_paths):
            raise ImportError("None of the files could be imported.")

        warning(f"Could not import {len(file_paths_missing_idx)} out of {len(file_paths)} files.")

        file_paths_missing = [file_paths[idx] for idx in file_paths_missing_idx]
        file_path_keep_boolean = [(file_path not in file_paths_missing) for file_path in file_paths]
        transcript_data = transcript_data.isel(
            file_path=file_path_keep_boolean,
            file=file_path_keep_boolean,
            drop=True,
        )

    if data_type == "stringtie" and read_length is not None:
        # we need to convert the coverage to counts
        # TODO: check if this works
        transcript_data["counts"] = transcript_data["counts"] * transcript_data["length"] / read_length

    if return_transcript_data:
        # return the transcript-level expression if requested
        if counts_from_abundance is not None:
            # convert the transcript counts to the desired count type
            transcript_data["counts"] = convert_abundance_to_counts(
                transcript_data["counts"],
                transcript_data["abundance"],
                transcript_data["length"],
                counts_from_abundance,
            )

        if ignore_after_bar:
            # change the transcript_id to only include the part before the bar for the coords
            transcript_data.coords["transcript_id"] = [
                transcript_id.split("|")[0] for transcript_id in transcript_data.coords["transcript_id"].values
            ]

        return transcript_data

    # convert to gene-level expression
    log(25, "Converting transcript-level expression to gene-level expression.")
    gene_expression: Union[xr.Dataset, ad.AnnData]
    gene_expression = convert_transcripts_to_genes(
        transcript_data,
        transcript_gene_map,
        ignore_after_bar=ignore_after_bar,
        ignore_transcript_version=ignore_transcript_version,
        counts_from_abundance=counts_from_abundance,
        biotype_filter=biotype_filter,
    )

    if output_type == "anndata":
        # create an AnnData object
        gene_expression = ad.AnnData(
            X=gene_expression["counts"].values.T,
            obs=pd.DataFrame(index=gene_expression.coords["file_path"].values),
            var=pd.DataFrame(index=gene_expression["gene_id"].values),
            obsm={
                "length": gene_expression["length"].values.T,
                "abundance": gene_expression["abundance"].values.T,
            },
        )

    if save_path is not None:
        if not isinstance(save_path, Path):
            save_path = Path(save_path)

        log(25, f"Saving the gene-level expression to: {save_path}.")

        if output_format == "h5ad":
            if not isinstance(gene_expression, ad.AnnData):
                raise ValueError("The output type must be AnnData to save as an h5ad file.")

            if save_path.suffix != ".h5ad":
                warning("The file extension of the `save_path` is not `.h5ad`. Appending the extension.")
                save_path = save_path.with_suffix(".h5ad")

            gene_expression.write(save_path)

        else:
            if isinstance(gene_expression, ad.AnnData):
                count_data = gene_expression.to_df().T
            else:
                count_data = gene_expression["counts"].to_pandas().values

            df_gene_data = pd.DataFrame(
                data=count_data,
                index=gene_expression["gene_id"],
                columns=gene_expression.coords["file_path"].values,
            )
            df_gene_data.sort_index(inplace=True)
            df_gene_data.to_csv(save_path, index=True, header=True, quoting=2)

    # end the timer
    end_time = time()
    log(25, f"Finished the import in {end_time - start_time:.2f} seconds.")

    if return_data:
        return gene_expression

    return None
