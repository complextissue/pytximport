from logging import log, warning
from pathlib import Path
from time import time
from typing import TYPE_CHECKING, Any, Callable, Dict, List, Literal, Optional, Union

import anndata as ad
import numpy as np
import pandas as pd
import xarray as xr
from tqdm import tqdm

from ..importers import read_kallisto, read_piscem, read_rsem, read_salmon, read_tsv
from ..utils import (
    convert_abundance_to_counts,
    convert_counts_to_tpm,
    convert_transcripts_to_genes,
    filter_by_biotype,
    get_median_length_over_isoform,
    remove_transcript_version,
    summarize_rsem_gene_data,
)

if TYPE_CHECKING:
    from summarizedexperiment import SummarizedExperiment


def tximport(
    file_paths: Union[List[str], List[Path]],
    data_type: Literal["kallisto", "salmon", "sailfish", "oarfish", "piscem", "stringtie", "rsem", "tsv"] = "salmon",
    transcript_gene_map: Optional[Union[pd.DataFrame, Union[str, Path]]] = None,
    counts_from_abundance: Optional[Literal["scaled_tpm", "length_scaled_tpm", "dtu_scaled_tpm"]] = None,
    gene_level: bool = False,
    return_transcript_data: bool = False,
    inferential_replicates: bool = False,
    inferential_replicate_transformer: Optional[Callable] = None,
    inferential_replicate_variance: bool = False,
    ignore_transcript_version: bool = True,
    ignore_after_bar: bool = True,
    id_column: Optional[str] = None,
    counts_column: Optional[str] = None,
    length_column: Optional[str] = None,
    abundance_column: Optional[str] = None,
    custom_importer: Optional[Callable] = None,
    existence_optional: bool = False,
    read_length: Optional[int] = None,
    # Arguments exclusive to the pytximport implementation
    output_type: Literal["xarray", "anndata", "summarizedexperiment"] = "anndata",
    output_format: Literal["csv", "h5ad", "summarizedexperiment"] = "csv",
    output_path: Optional[Union[str, Path]] = None,
    output_path_overwrite: bool = False,
    return_data: bool = True,
    biotype_filter: Optional[List[str]] = None,
) -> Union[xr.Dataset, ad.AnnData, "SummarizedExperiment", None]:
    """Import transcript-level quantification files and convert them to gene-level expression estimates.

    Basic usage:

    .. code-block:: python

        from pytximport import tximport

        txi = tximport(
            ["quant_1.sf", "quant_2.sf"],
            data_type="salmon",
            transcript_gene_map="transcript_to_gene_map.tsv",
            counts_from_abundance="length_scaled_tpm",
        )

    Args:
        file_paths (List[Union[str, Path]]): The paths to the quantification files.
        data_type (Literal["kallisto", "salmon", "sailfish", "oarfish", "piscem", "stringtie", "rsem", "tsv"]): The type
            of quantification files. Defaults to "salmon".
        transcript_gene_map (Optional[Union[pd.DataFrame, Union[str, Path]], optional): The mapping from transcripts to
            genes. Has to contain two columns: `transcript_id` and `gene_id`. If you provide a path to a file, it has to
            be either a tab-separated (.tsv) or comma-separated (.csv) file with a header. Defaults to None.
        counts_from_abundance (Optional[Literal["scaled_tpm", "length_scaled_tpm", "dtu_scaled_tpm"]], optional):
            Whether to calculate count estimates based on the abundance. When using scaled_tpm or length_scaled_tpm the
            counts no longer correlate with the the average transcript length per sample. In those cases, the length
            offset matrix should not be used for downstream analysis. Note, that this does not normalize the sequencing
            depth, only the difference in transcript length. When using the gene-summarized counts and not count
            estimates based on the abundance, the length offset matrix included in the output from this function should
            be used for downstream analysis. If your downstream analysis tool does not support the length offset matrix,
            you should probably use `length_scaled_tpm` for gene-level analysis. For transcript-level analysis, we
            recommend that you use `scaled_tpm` or `dtu_scaled_tpm`. For further guidance on transcript-level analysis,
            please refer to: https://doi.org/10.12688/f1000research.15398.3.
            Defaults to None.
        gene_level (bool, optional): Whether the input files are at the gene level. This is only the case for some RSEM
            quantification files. Defaults to False.
        return_transcript_data (bool, optional): Whether to return the transcript-level expression. Defaults to False.
        inferential_replicates (bool, optional): Whether to parse and include inferential replicates in the output.
            If you want to recalculate the counts from inferential replicates, please set this option to True and
            provide a `inferential_replicate_transformer`. Defaults to False.
        inferential_replicate_transformer (Optional[Callable], optional): A custom function to transform the inferential
            replicates. Defaults to None.
        inferential_replicate_variance (bool, optional): Whether to return the variance of the inferential replicates.
            Defaults to False.
        ignore_transcript_version (bool, optional): Whether to ignore the transcript version. Defaults to True.
        ignore_after_bar (bool, optional): Whether to split the transcript id after the bar character (`|`).
            Defaults to True.
        id_column (Optional[str], optional): The column name for the transcript id. Defaults to None.
        counts_column (Optional[str], optional): The column name for the counts. Defaults to None.
        length_column (Optional[str], optional): The column name for the length. Defaults to None.
        abundance_column (Optional[str], optional): The column name for the abundance. Defaults to None.
        custom_importer (Optional[Callable], optional): A custom importer function. Defaults to None.
        existence_optional (bool, optional): Whether the existence of the files is optional. Defaults to False.
        read_length (Optional[int], optional): The read length for the stringtie quantification. Defaults to None.
        output_type (Literal["xarray", "anndata"], optional): The type of output. Defaults to "anndata".
        output_format (Literal["csv", "h5ad"], optional): The type of output file. Defaults to "csv".
        output_path (Optional[Union[str, Path]], optional): The path to save the gene-level expression.
            Defaults to None.
        output_path_overwrite (bool, optional): Whether to overwrite the save path if it already exists.
            Defaults to False.
        return_data (bool, optional): Whether to return the gene-level expression. Defaults to True.
        biotype_filter (List[str], optional): Filter the transcripts by biotype, including only those provided. Enables
            post-hoc filtering of the data based on the biotype of the transcripts. Assumes that the biotype is present
            in the transcript_id of the data, bar-separated. Defaults to None.

    Returns:
        Union[xr.Dataset, ad.AnnData, SummarizedExperiment, None]: The estimated gene-level or transcript-level
            expression data if `return_data` is True, else None.
    """
    # Start a timer
    log(25, "Starting the import.")
    start_time = time()

    # Needed for faster groupby operations
    xr.set_options(use_flox=True)

    assert data_type in [
        "kallisto",
        "salmon",
        "sailfish",
        "oarfish",
        "piscem",
        "stringtie",
        "rsem",
        "tsv",
    ], "Only kallisto, salmon/sailfish, oarfish, piscem, stringtie, RSEM and tsv quantification files are supported."

    # Read the transcript to gene mapping
    if isinstance(transcript_gene_map, str) or isinstance(transcript_gene_map, Path):
        transcript_gene_map = Path(transcript_gene_map)
        if not transcript_gene_map.exists():
            raise FileNotFoundError(f"The transcript to gene mapping does not exist: {transcript_gene_map}")

        try:
            transcript_gene_map = pd.read_csv(
                transcript_gene_map,
                header=0,
                engine="c",
                sep=("," if transcript_gene_map.suffix == ".csv" else "\t"),
                usecols=["transcript_id", "gene_id"],
                dtype={"transcript_id": str, "gene_id": str},
            )
        except Exception as exception:
            raise ValueError(f"Could not read the transcript to gene mapping: {exception}")

    if transcript_gene_map is not None:
        # Assert that transcript_id and gene_id are present in the mapping
        assert isinstance(transcript_gene_map, pd.DataFrame), "The mapping must be a DataFrame."
        assert "transcript_id" in transcript_gene_map.columns, "The mapping does not contain a `transcript_id` column."
        assert "gene_id" in transcript_gene_map.columns, "The mapping does not contain a `gene_id` column."

        # Check whether the mapping contains duplicates
        if transcript_gene_map.duplicated(subset=["transcript_id", "gene_id"]).any():
            warning("The transcript to gene mapping contains duplicates. Removing duplicates.")
            transcript_gene_map = transcript_gene_map.drop_duplicates(
                subset=["transcript_id", "gene_id"],
                keep="first",
            )

    # Assert that return_transcript_data is True if transcript_gene_map is None
    if transcript_gene_map is None:
        assert (
            return_transcript_data or gene_level
        ), "A transcript to gene mapping must be provided when summarizing transcripts to genes."

    if gene_level and data_type != "rsem":
        raise ValueError("Gene-level imports are only available for RSEM quantification files.")

    # Match the data type with the required importer
    importer: Callable
    importer_kwargs: Dict[str, Any] = {}

    if data_type == "kallisto":
        importer = read_kallisto
    elif data_type == "salmon" or data_type == "sailfish":
        importer = read_salmon
    elif data_type == "rsem":
        if gene_level and return_transcript_data:
            raise ValueError("Cannot return transcript-level data from gene-level counts")

        if gene_level and counts_from_abundance is not None:
            raise ValueError("Count estimation from abundances requires transcript-level data.")

        id_column = ("gene_id" if gene_level else "transcript_id") if id_column is None else id_column

        importer_kwargs["gene_level"] = gene_level
        importer = read_rsem
    elif data_type == "oarfish":
        # Use the default column names if not provided
        id_column = "tname" if id_column is None else id_column
        counts_column = "num_reads" if counts_column is None else counts_column
        length_column = "len" if length_column is None else length_column
        abundance_column = "num_reads" if abundance_column is None else abundance_column

        importer = read_piscem
    elif data_type == "piscem":
        warning(
            (
                "Assuming a piscem-infer .quant file with columns: target_name, ecount, eeln, tpm. "
                "This differs from the assumed columns in the original tximport package. "
                "If you encounter issues, please provide the column names explicitly."
            )
        )

        # Column names based on https://piscem-infer.readthedocs.io/en/latest/formats.html
        id_column = "target_name" if id_column is None else id_column
        counts_column = "ecount" if counts_column is None else counts_column
        length_column = "eeln" if length_column is None else length_column
        abundance_column = "tpm" if abundance_column is None else abundance_column

        importer = read_piscem
    elif data_type == "stringtie":
        if read_length is None:
            raise ValueError("The read_length must be provided for stringtie quantification files.")

        if any(Path(file_path).suffix == ".gtf" for file_path in file_paths):
            raise ValueError("The input file type is a GTF file. Please provide a tab-separated file instead.")

        id_column = "t_name" if id_column is None else id_column
        counts_column = "cov" if counts_column is None else counts_column
        length_column = "length" if length_column is None else length_column
        abundance_column = "FPKM" if abundance_column is None else abundance_column

        importer = read_tsv
    elif data_type == "tsv":
        # Check that all of id_column counts_column length_column abundance_column are provided
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

    # Overwrite the importer if a custom importer is provided
    if custom_importer is not None:
        importer = custom_importer

    # Create the importer kwargs based on the provided column names
    importer_kwargs.update(
        {
            "id_column": id_column,
            "counts_column": counts_column,
            "length_column": length_column,
            "abundance_column": abundance_column,
        }
    )

    # List is required to create a copy so that the original dictionary can be changed
    for key, value in list(importer_kwargs.items()):
        if value is None:
            del importer_kwargs[key]

    if data_type in ["salmon", "sailfish", "kallisto", "piscem", "oarfish"] and inferential_replicates:
        # Add information about the inferential replicates to the importer kwargs
        if data_type == "salmon":
            importer_kwargs["aux_dir_name"] = "aux_info"
        elif data_type == "sailfish":
            # TODO: this may break in newer versions of sailfish (>0.10.1), when no explicit auxDir is provided
            importer_kwargs["aux_dir_name"] = "aux"

        # Set whether the counts and abundances will be recomputed from the inferential replicates
        recompute_counts = inferential_replicate_transformer is not None

        # Pass whether to load inferential replicates to the importer
        importer_kwargs["inferential_replicates"] = inferential_replicates
        importer_kwargs["recompute_counts"] = recompute_counts
    elif inferential_replicates:
        warning("Inferential replicates are not supported for this data type.")
        inferential_replicates = False
    else:
        recompute_counts = False

    transcript_data: Optional[xr.Dataset] = None
    file_paths_missing_idx: List[int] = []

    # Iterate through the files, unless they are RSEM gene level counts
    if not (gene_level and data_type == "rsem"):
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

            # The check for transcript_data is needed in case the import of the first file fails
            if file_idx == 0 or transcript_data is None:
                empty_array = np.zeros((len(transcript_data_sample["transcript_id"]), len(file_paths)))

                abundance = xr.DataArray(data=empty_array.copy(), dims=["transcript_id", "file"])
                counts = xr.DataArray(data=empty_array.copy(), dims=["transcript_id", "file"])
                length = xr.DataArray(data=empty_array.copy(), dims=["transcript_id", "file"])

                data_vars = {
                    "abundance": abundance,
                    "counts": counts,
                    "length": length,
                }

                if inferential_replicates:
                    if transcript_data_sample["inferential_replicates"] is None:
                        raise ValueError(
                            f"The quantification file does not contain inferential replicates: {file_path}"
                        )

                    bootstrap_count = transcript_data_sample["inferential_replicates"]["replicates"].shape[1]
                    data_vars["inferential_replicates"] = xr.DataArray(
                        data=np.zeros((len(transcript_data_sample["transcript_id"]), bootstrap_count, len(file_paths))),
                        dims=["transcript_id", "bootstrap", "file"],
                    )

                    if inferential_replicate_variance:
                        data_vars["variance"] = xr.DataArray(data=empty_array.copy(), dims=["transcript_id", "file"])

                transcript_data = xr.Dataset(
                    data_vars,
                    coords={
                        "transcript_id": transcript_data_sample["transcript_id"],
                        "file_path": list(file_paths),
                    },
                )

            # Add the transcript-level expression to the array
            # Counts and abundance depend on whether inferential replicates are used and are added later
            transcript_data["length"].loc[{"file": file_idx}] = transcript_data_sample["length"]

            if inferential_replicates:
                if transcript_data_sample["inferential_replicates"] is None:
                    raise ValueError(f"The quantification file does not contain inferential replicates: {file_path}")

                transcript_data["inferential_replicates"].loc[{"file": file_idx}] = transcript_data_sample[
                    "inferential_replicates"
                ]["replicates"]

                if inferential_replicate_variance:
                    transcript_data["variance"].loc[{"file": file_idx}] = transcript_data_sample[
                        "inferential_replicates"
                    ]["variance"]

            if recompute_counts:
                # Use the provided function to recompute the counts based on the inferential replicates
                transcript_data["counts"].loc[{"file": file_idx}] = inferential_replicate_transformer(  # type: ignore
                    transcript_data_sample["inferential_replicates"]["replicates"],
                )

                # Recompute the abundance
                transcript_data["abundance"].loc[{"file": file_idx}] = convert_counts_to_tpm(
                    transcript_data["counts"].loc[{"file": file_idx}].data,
                    transcript_data["length"].loc[{"file": file_idx}].data,
                )
            else:
                # We only need to add the abundance and counts if we are not recalculating them
                transcript_data["abundance"].loc[{"file": file_idx}] = transcript_data_sample["abundance"]
                transcript_data["counts"].loc[{"file": file_idx}] = transcript_data_sample["counts"]
    else:
        transcript_data, file_paths_missing_idx = summarize_rsem_gene_data(
            file_paths,
            importer,
            importer_kwargs,
            existence_optional,
        )

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
        # We need to convert the coverage to counts
        # TODO: add a unit test that checks for agreement with tximport
        transcript_data["counts"] = transcript_data["counts"] * transcript_data["length"] / read_length

    if biotype_filter is not None:
        transcript_data = filter_by_biotype(
            transcript_data, biotype_filter, id_column=("gene_id" if gene_level else "transcript_id")
        )

    # Remove appended gene names after underscore for RSEM data for both transcript and gene ids
    if (
        data_type == "rsem"
        and ignore_after_bar
        and (
            (gene_level and transcript_data.coords["gene_id"].values[0].count("_") > 0)
            or (not gene_level and transcript_data.coords["transcript_id"].values[0].count("_") > 0)
        )
    ):
        warning(
            (
                "RSEM seems to have been run with `--append-names`. "
                "Removing the appended names. "
                "Set `ignore_after_bar` to False to keep the appended names."
            )
        )
        if not gene_level:
            transcript_data.coords["transcript_id"] = [
                transcript_id.split("_")[0] for transcript_id in transcript_data.coords["transcript_id"].values
            ]
        else:
            transcript_data.coords["gene_id"] = [
                gene_id.split("_")[0] for gene_id in transcript_data.coords["gene_id"].values
            ]

    result: Union[xr.Dataset, ad.AnnData]
    result_index = "transcript_id" if return_transcript_data else "gene_id"
    data_id_column = "transcript_id" if not gene_level else "gene_id"

    if ignore_after_bar:
        # Ignore the part of the transcript ID after the bar
        transcript_data.coords[data_id_column] = [
            id.split("|")[0] for id in transcript_data.coords[data_id_column].values
        ]

    if ignore_transcript_version:
        transcript_ids = transcript_data.coords["transcript_id"].values if not gene_level else None

        # Ignore the transcript version in both the data and the transcript gene map
        transcript_data, transcript_gene_map, transcript_ids = remove_transcript_version(  # type: ignore
            transcript_data,
            transcript_gene_map,
            transcript_ids,  # type: ignore
            id_column=data_id_column,
        )

    if return_transcript_data:
        if counts_from_abundance is not None:
            if counts_from_abundance == "length_scaled_tpm":
                raise ValueError("The `length_scaled_tpm` option is not supported for transcript-level expression.")

            if counts_from_abundance == "dtu_scaled_tpm":
                if transcript_gene_map is None:
                    raise ValueError("A transcript to gene mapping must be provided for `dtu_scaled_tpm`.")

                log(25, "Calculating median gene length over isoforms.")
                transcript_data = get_median_length_over_isoform(
                    transcript_data,
                    transcript_gene_map,
                )
                counts_from_abundance = "length_scaled_tpm"
                length_key = "median_isoform_length"
            else:
                length_key = "length"

            # Convert the transcript counts to the desired count type
            log(25, "Recreating transcript counts from abundances.")
            transcript_data["counts"] = convert_abundance_to_counts(
                transcript_data["counts"],
                transcript_data["abundance"],
                transcript_data[length_key],
                counts_from_abundance,
            )

        result = transcript_data
    elif gene_level:
        if counts_from_abundance is not None:
            raise ValueError("Count estimation from abundances requires transcript-level data.")

        result = transcript_data
    else:
        # Convert to gene-level expression
        log(25, "Converting transcript-level expression to gene-level expression.")
        if counts_from_abundance == "dtu_scaled_tpm":
            raise ValueError("The `dtu_scaled_tpm` option is not supported for gene-level expression.")

        result = convert_transcripts_to_genes(
            transcript_data,
            transcript_gene_map,  # type: ignore
            counts_from_abundance=counts_from_abundance,
        )
        result_index = "gene_id"

    if output_path is not None:
        if not isinstance(output_path, Path):
            output_path = Path(output_path)

        if output_path.suffix == ".csv" and output_format == "h5ad":
            warning(
                "The file extension of the `output_path` is `.csv` but the output format is `.h5ad`. "
                "Changing the output format to `.csv`."
            )
            output_format = "csv"

        if output_format == "h5ad" and output_type != "anndata":
            warning(
                "The output format is h5ad but the output type is not anndata. Changing the output type to anndata."
            )
            output_type = "anndata"

        if output_format == "summarizedexperiment" and output_type != "summarizedexperiment":
            warning(
                "The output format is summarizedexperiment but the output type is not summarizedexperiment. "
                "Changing the output type to summarizedexperiment."
            )
            output_type = "summarizedexperiment"

    if output_type == "anndata":
        obsm = {
            "length": result["length"].values.T,
            "abundance": result["abundance"].values.T,
        }
        uns: Dict[str, Any] = {
            "counts_from_abundance": counts_from_abundance,
        }

        if inferential_replicates:
            if "variance" in result.data_vars:
                obsm["variance"] = result["variance"].values.T

            uns["inferential_replicates"] = np.moveaxis(result["inferential_replicates"].values, -1, 0)

        result = ad.AnnData(
            X=result["counts"].values.T,
            obs=pd.DataFrame(index=result.coords["file_path"].values),
            var=pd.DataFrame(index=result[result_index].values),
            obsm=obsm,
            uns=uns,
        )
    elif output_type == "summarizedexperiment":
        try:
            from biocframe import BiocFrame
            from summarizedexperiment import SummarizedExperiment
        except ImportError:
            raise ImportError(
                "To export data as SummarizedExperiment, please install the optional dependencies from BiocPy. "
                "E.g.: `pip install pytximport[biocpy]`"
            )

        warning("Support for the SummarizedExperiment output type is experimental.")

        summarized_experiment_metadata = {
            "counts_from_abundance": counts_from_abundance,
            "length": result["length"].values,
            "abundance": result["abundance"].values,
        }

        if inferential_replicates:
            if "variance" in result.data_vars:
                summarized_experiment_metadata["variance"] = result["variance"].values

            summarized_experiment_metadata["inferential_replicates"] = result["inferential_replicates"].values

        result = SummarizedExperiment(
            assays={"counts": result["counts"].values},
            row_data=BiocFrame(row_names=result[result_index].values),
            column_data=BiocFrame(row_names=result.coords["file_path"].values),
            metadata=summarized_experiment_metadata,
        )

    if output_path is not None:
        if output_path.exists() and not output_path_overwrite:
            raise FileExistsError(
                f"The file already exists: {output_path}. Set `output_path_overwrite` to True to overwrite."
            )

        if not output_path.parent.exists():
            output_path.parent.mkdir(parents=True)

        log(25, f"Saving the gene-level expression to: {output_path}.")

        if output_format == "h5ad":
            if not isinstance(result, ad.AnnData):
                raise ValueError("The output type must be AnnData to save as an h5ad file.")

            if output_path.suffix != ".h5ad":
                warning("The file extension of the `output_path` is not `.h5ad`. Appending the extension.")
                output_path = output_path.with_suffix(".h5ad")

            result.write(output_path)
        elif output_format == "summarizedexperiment":
            try:
                from dolomite_base import save_object
            except Exception as e:
                raise ImportError(
                    "Please install the optional SummarizedExperiment dependencies from BiocPy.\n"
                    "`pip install pytximport[biocpy] \n\n"
                ) from e

            if not isinstance(result, SummarizedExperiment):
                raise ValueError("The output type must be 'summarizedexperiment' to save as file.")

            save_object(result, str(output_path))
        else:
            if output_type == "summarizedexperiment":
                # to avoid a top level import to the summarizedexperiment package.
                try:
                    count_data = result.assay("counts")
                except Exception:
                    raise Exception(
                        "Could not convert the SummarizedExperiment object to a DataFrame, to save as a CSV. "
                        "Please choose `summarizedexperiment` as the output format instead."
                    )

                df_gene_data = pd.DataFrame(
                    data=count_data,
                    index=result.get_row_names(),
                    columns=result.get_column_names(),
                )
                df_gene_data.sort_index(inplace=True)
                df_gene_data.to_csv(output_path, index=True, header=True, quoting=2)
            else:
                if isinstance(result, ad.AnnData):
                    try:
                        count_data = result.to_df().T
                    except Exception:
                        raise Exception(
                            "Could not convert the AnnData object to a DataFrame, to save as a CSV. "
                            "Please choose `.h5ad` as the output format instead or `xarray` as the output type."
                        )
                else:
                    count_data = result["counts"].to_pandas().values

                df_gene_data = pd.DataFrame(
                    data=count_data,
                    index=(result[result_index] if output_type != "anndata" else result.var.index),
                    columns=(result.coords["file_path"].values if output_type != "anndata" else result.obs.index),
                )
                df_gene_data.sort_index(inplace=True)
                df_gene_data.to_csv(output_path, index=True, header=True, quoting=2)

    # End the timer
    log(25, f"Finished the import in {time() - start_time:.2f} seconds.")

    if return_data:
        return result

    return None
