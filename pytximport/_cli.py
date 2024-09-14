"""Expose the tximport function as a command-line tool."""

from logging import basicConfig

import click
import numpy as np

from .core import tximport


@click.command()
@click.option(
    "-i",
    "--file_paths",
    "--file-paths",
    type=click.Path(exists=False),
    multiple=True,
    help="The path to an quantification file. To provide multiple input files, use `-i input1.sf -i input2.sf ...`.",
    required=True,
)
@click.option(
    "-t",
    "--data_type",
    "--data-type",
    type=click.Choice(["kallisto", "salmon", "sailfish", "oarfish", "piscem", "stringtie", "rsem", "tsv"]),
    help="The type of quantification files.",
    required=True,
)
@click.option(
    "-m",
    "--transcript_gene_map",
    "--transcript-gene-map",
    type=click.Path(exists=True),
    help=(
        "The path to the transcript to gene map. Either a tab-separated (.tsv) or comma-separated (.csv) file. "
        "Expected column names are `transcript_id` and `gene_id`."
    ),
)
@click.option(
    "-c",
    "--counts_from_abundance",
    "--counts-from-abundance",
    type=click.Choice(["scaled_tpm", "length_scaled_tpm", "dtu_scaled_tpm"]),
    help=(
        "The method to calculate the counts from the abundance. Leave empty to use counts. "
        "For differential gene expression analysis, we recommend using `length_scaled_tpm`. "
        "For differential transcript expression analysis, we recommend using `scaled_tpm`. "
        "For differential isoform usage analysis, we recommend using `dtu_scaled_tpm`."
    ),
)
@click.option(
    "-o",
    "--output_path",
    "--save-path",
    type=click.Path(),
    help="The output path to save the resulting counts to.",
    required=True,
)
@click.option(
    "-of",
    "--output_format",
    "--output-format",
    type=click.Choice(["csv", "h5ad"]),
    help="The format of the output file.",
)
@click.option(
    "-ow",
    "--output_path_overwrite",
    "--save-path-overwrite",
    is_flag=True,
    help="Provide this flag to overwrite an existing file at the output path.",
)
@click.option(
    "--ignore_after_bar",
    "--ignore-after-bar",
    type=bool,
    default=True,
    help="Whether to split the transcript id after the bar character (`|`).",
)
@click.option(
    "--ignore_transcript_version",
    "--ignore-transcript-version",
    type=bool,
    default=True,
    help="Whether to ignore the transcript version.",
)
@click.option(
    "-gl",
    "--gene_level",
    "--gene-level",
    is_flag=True,
    help="Provide this flag when importing gene-level counts from RSEM files.",
)
@click.option(
    "-tx",
    "--return_transcript_data",
    "--return-transcript-data",
    is_flag=True,
    help=(
        "Provide this flag to return transcript-level instead of gene-summarized data. "
        "Incompatible with gene-level input and `counts_from_abundance=length_scaled_tpm`."
    ),
)
@click.option(
    "-ir" "--inferential_replicates",
    "--inferential-replicates",
    is_flag=True,
    help="Provide this flag to make use of inferential replicates. Will use the median of the inferential replicates.",
)
@click.option(
    "-id",
    "--id_column",
    "--id-column",
    type=str,
    help="The column name for the transcript id.",
)
@click.option(
    "-counts",
    "--counts_column",
    "--counts-column",
    type=str,
    help="The column name for the counts.",
)
@click.option(
    "-length",
    "--length_column",
    "--length-column",
    type=str,
    help="The column name for the length.",
)
@click.option(
    "-tpm",
    "--abundance_column",
    "--abundance-column",
    type=str,
    help="The column name for the abundance.",
)
@click.option(
    "--existence_optional",
    "--existence-optional",
    is_flag=True,
    help="Whether the existence of the files is optional.",
)
def cli(  # type: ignore
    **kwargs,
) -> None:
    """Call the tximport function via the command line.

    You can view the available options by running `pytximport --help`.

    .. code-block:: bash

        pytximport --help

    For detailed information on pytximport's functionality, please refer to the README and online documentation.

    Args:
        **kwargs: The keyword arguments to pass to the tximport function.

    Returns:
        None
    """
    # Add return_data to the kwargs with a default value of False
    kwargs["return_data"] = False
    kwargs["output_type"] = "anndata"
    kwargs["inferential_replicate_transformer"] = lambda x: np.median(x, axis=1)

    # Set the logging level
    basicConfig(level=25, format="%(asctime)s: %(message)s")

    tximport(**kwargs)  # type: ignore
