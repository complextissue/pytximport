"""Collections of biotype filters for use with the `biotype_filter` option of `pytximport`.

You can also pass your own list of biotypes to filter by. For example:

.. code-block:: python

    from pytximport import tximport

    txi = tximport(
        file_paths,
        data_type,
        transcript_gene_map,
        counts_from_abundance,
        biotype_filter=["protein_coding", "nonsense_mediated_decay"],
    )

For a list of all biotypes, refer to GENCODE's biotype documentation: https://www.gencodegenes.org/pages/biotypes.html
If you use a different reference/annotation, available biotypes may differ.
"""

from typing import List

GENCODE_PROTEIN_CODING: List[str] = [
    "protein_coding",
    "nonsense_mediated_decay",
    "non_stop_decay",
    "IG_C_gene",
    "IG_D_gene",
    "IG_J_gene",
    "IG_LV_gene",
    "IG_V_gene",
    "TR_C_gene",
    "TR_J_gene",
    "TR_V_gene",
    "TR_D_gene",
    "polymorphic_pseudogene",
    "protein_coding_LoF",
]
