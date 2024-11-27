"""Test the cli function."""

from pathlib import Path
from typing import List

import anndata as ad
import numpy as np
import xarray as xr

from pytximport import tximport
from pytximport.utils import create_transcript_gene_map, filter_by_biotype


def test_biotype_filter(
    salmon_multiple_files: List[Path],
) -> None:
    """Test filtering an existing transcript- or gene-level AnnData object or xaarray Dataset by biotype."""
    transcript_gene_map_mouse = create_transcript_gene_map(
        species="mouse",
        target_field=["external_gene_name", "gene_biotype"],
    )

    assert "gene_biotype" in transcript_gene_map_mouse.columns, "The gene_biotype column is missing."

    for transcript_level in [True, False]:
        for output_type in ["anndata", "xarray"]:
            result = tximport(
                salmon_multiple_files,
                "salmon",
                transcript_gene_map_mouse,
                output_type=output_type,
                return_transcript_data=transcript_level,
            )

            for recalculate_abundance in [True, False]:
                id_column = "transcript_id" if transcript_level else "gene_id"

                result_filtered = filter_by_biotype(
                    result,
                    transcript_gene_map_mouse,
                    biotype_filter=["protein_coding"],
                    id_column=id_column,
                    recalculate_abundance=recalculate_abundance,
                )

                if output_type == "anndata":
                    assert isinstance(result_filtered, ad.AnnData)
                    assert len(result_filtered.var_names) > 0, "No genes were retained."
                    assert len(result_filtered.var_names) < len(result.var_names), "All genes were retained."

                    if recalculate_abundance:
                        np.testing.assert_allclose(
                            result_filtered.obsm["abundance"].sum(axis=1),
                            result.obsm["abundance"].sum(axis=1),
                            rtol=1e-4,
                            atol=1,
                        )

                else:
                    assert isinstance(result_filtered, xr.Dataset)
                    assert len(result_filtered.coords[id_column]) > 0, "No genes were retained."
                    assert len(result_filtered.coords[id_column]) < len(
                        result.coords[id_column]
                    ), "All genes were retained."

                    if recalculate_abundance:
                        np.testing.assert_allclose(
                            result_filtered["abundance"].sum(axis=0),
                            result["abundance"].sum(axis=0),
                            rtol=1e-4,
                            atol=1,
                        )
