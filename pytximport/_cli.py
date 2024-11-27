"""Expose the tximport function as a command-line tool."""

from logging import basicConfig, log, warning
from pathlib import Path

import click
import numpy as np
from click_default_group import DefaultGroup

from .core import tximport
from .utils import create_transcript_gene_map_from_annotation


@click.group(
    cls=DefaultGroup,
    default="run",
    default_if_no_args=True,
    help="Welcome to the pytximport command-line interface for importing transcript-level quantification files.",
)
@click.pass_context
def cli(  # type: ignore  # pragma: no cover
    ctx: click.Context,
):
    """Welcome to the pytximport command-line interface for importing transcript-level quantification files."""
    pass


@cli.command(
    no_args_is_help=True,
)
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
def run(  # type: ignore  # pragma: no cover
    **kwargs,
) -> None:
    """Call the tximport function via the command line."""
    basicConfig(level=25, format="%(asctime)s: %(message)s")

    # Add return_data to the kwargs with a default value of False
    kwargs["return_data"] = False
    kwargs["output_type"] = "anndata"
    kwargs["inferential_replicate_transformer"] = lambda x: np.median(x, axis=1)

    tximport(**kwargs)  # type: ignore


@cli.command(
    no_args_is_help=True,
)
@click.option(
    "-i",
    "--input_file",
    "--input",
    type=click.Path(exists=True),
    help="The path to the annotation GTF file.",
    required=True,
)
@click.option(
    "-o",
    "--output_file",
    "--output",
    type=click.Path(),
    help="The output path to save the resulting transcript-to-gene mapping file to.",
    required=True,
)
@click.option(
    "-ow",
    "--output_path_overwrite",
    "--save-path-overwrite",
    is_flag=True,
    help="Provide this flag to overwrite an existing file at the output path.",
)
@click.option(
    "--source-field",
    "--source_field",
    type=str,
    help="The annotation field to use as the source in the mapping file.",
    required=False,
)
@click.option(
    "--target-field",
    "--target_field",
    type=str,
    multiple=True,
    help="The annotation field(s) to use as the target in the mapping file.",
    required=False,
)
@click.option(
    "--keep-biotype",
    "--keep_biotype",
    is_flag=True,
    help="Provide this flag to keep the gene_biotype column as an additional column in the mapping file.",
)
def create_map(  # type: ignore  # pragma: no cover
    **kwargs,
) -> None:
    """Create a transcript-to-gene mapping file via the command line."""
    basicConfig(level=25, format="%(asctime)s: %(message)s")
    log(25, "Creating a transcript-to-gene mapping file.")

    if isinstance(kwargs["target_field"], tuple):
        kwargs["target_field"] = list(kwargs["target_field"])

    df = create_transcript_gene_map_from_annotation(
        kwargs["input_file"],
        source_field=kwargs["source_field"] if kwargs["source_field"] else "transcript_id",
        target_field=kwargs["target_field"] if kwargs["target_field"] else "gene_id",
        keep_biotype=kwargs["keep_biotype"],
    )
    log(25, "Created the transcript-to-gene mapping file. Saving the file...")

    output_file = Path(kwargs["output_file"])
    if not output_file.exists() or kwargs["output_path_overwrite"]:
        df.to_csv(
            kwargs["output_file"],
            sep=("," if kwargs["output_file"].endswith(".csv") else "\t"),
            index=False,
        )
        log(25, f"Saved the transcript-to-gene mapping file to {kwargs['output_file']}.")
    else:
        warning(
            f"Could not save the transcript-to-gene mapping file. File already exists at {kwargs['output_file']}. "
            "Use the `-ow` flag to overwrite."
        )
