"""Expose the tximport function as a command-line tool."""

from logging import basicConfig

import click

from .core import tximport


@click.command()
@click.option(
    "-i",
    "--file_paths",
    "--file-paths",
    type=click.Path(exists=False),
    multiple=True,
    help="The paths to the quantification files.",
    required=True,
)
@click.option(
    "-t",
    "--data_type",
    "--data-type",
    type=click.Choice(["kallisto", "salmon", "sailfish", "oarfish", "piscem", "stringtie", "tsv"]),
    help="The type of quantification file.",
    required=True,
)
@click.option(
    "-m",
    "--transcript_gene_map",
    "--transcript-gene-map",
    type=click.Path(exists=True),
    help="The path to the transcript to gene mapping file.",
)
@click.option(
    "-c",
    "--counts_from_abundance",
    "--counts-from-abundance",
    type=click.Choice(["scaled_tpm", "length_scaled_tpm", "dtu_scaled_tpm"]),
    help="The type of counts to convert to.",
)
@click.option(
    "-o",
    "--save_path",
    "--save-path",
    type=click.Path(),
    help="The path to save the gene-level expression.",
    required=True,
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
    "-tx",
    "--return_transcript_data",
    "--return-transcript-data",
    is_flag=True,
    help="Whether to return transcript-level instead of gene-summarized data.",
)
@click.option(
    "-ir" "--inferential_replicates",
    "--inferential-replicates",
    is_flag=True,
    help="Whether to make use of inferential replicates.",
)
@click.option(
    "-irv",
    "-inferential_replicate_variance",
    "--inferential-replicate-variance",
    is_flag=True,
    help="Whether to calculate the variance from the inferential replicates.",
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
    is_flag=True,
    help="Whether the existence of the files is optional.",
)
@click.option(
    "--output_type",
    "--output-type",
    type=click.Choice(["xarray", "anndata"]),
    help="The type of output to return.",
)
@click.option(
    "--output_format",
    "--output-format",
    type=click.Choice(["csv", "h5ad"]),
    help="The format of the output file.",
)
def cli(  # type: ignore
    **kwargs,
) -> None:
    """Convert transcript-level expression to gene-level expression."""
    # add return_data to the kwargs with a default value of False
    kwargs["return_data"] = False

    # set the logging level
    basicConfig(level=25, format="%(asctime)s: %(message)s")

    tximport(**kwargs)  # type: ignore
