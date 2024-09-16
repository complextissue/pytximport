"""Benchmark pytximport."""

from logging import INFO, getLogger

import numpy as np
import pandas as pd
import pyperf

from pytximport import tximport

# Load transcript-gene mapping (this can be done once outside the benchmark)
transcript_gene_mapping_human = pd.read_table(
    "../data/fabry_disease/transcript_gene_mapping_human.tsv",
    header=0,
    sep="\t",
)

# Define the files list
files = [
    "../data/fabry_disease/SRR16504309_wt/quant.sf",
    "../data/fabry_disease/SRR16504310_wt/quant.sf",
    "../data/fabry_disease/SRR16504311_ko/quant.sf",
    "../data/fabry_disease/SRR16504312_ko/quant.sf",
]


# Function to benchmark
def tximport_benchmark():
    """Benchmark pytximport.

    Returns:
        ad.AnnData: The AnnData object.
    """
    txi = tximport(
        files,
        "salmon",
        transcript_gene_mapping_human,
        inferential_replicates=True,
        inferential_replicate_transformer=lambda x: np.median(x, axis=1),
    )
    return txi


# Run the pyperf benchmark
if __name__ == "__main__":
    # Set log level to 25
    getLogger().setLevel(25)

    runner = pyperf.Runner()
    runner.bench_func("pytximport", tximport_benchmark)
