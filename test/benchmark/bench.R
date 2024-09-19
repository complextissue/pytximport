options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("tximport", quietly = TRUE)) {
    install.packages("BiocManager")
    BiocManager::install("tximport")
}

if (!requireNamespace("bench", quietly = TRUE)) {
    install.packages("bench")
}

if (!requireNamespace("readr", quietly = TRUE)) {
    install.packages("readr")
}

if (!requireNamespace("matrixStats", quietly = TRUE)) {
    install.packages("matrixStats")
}

# Disable multithreading (should only use a single thread anyway)
Sys.setenv("OMP_NUM_THREADS" = "1")
Sys.setenv("OPENBLAS_NUM_THREADS" = "1")
Sys.setenv("MKL_NUM_THREADS" = "1")
Sys.setenv("VECLIB_MAXIMUM_THREADS" = "1")
Sys.setenv("NUMEXPR_NUM_THREADS" = "1")

# Load libraries
library(tximport)
library(bench)
library(readr)
library(matrixStats)

# Define the file paths and transcript-gene mapping
transcript_gene_mapping <- read_tsv("../data/fabry_disease/transcript_gene_mapping_human.tsv")

files <- c(
    "../data/fabry_disease/SRR16504309_wt/quant.sf",
    "../data/fabry_disease/SRR16504310_wt/quant.sf",
    "../data/fabry_disease/SRR16504311_ko/quant.sf",
    "../data/fabry_disease/SRR16504312_ko/quant.sf"
)

# Create a function to benchmark
tximport_benchmark <- function() {
    txi <- tximport(
        files,
        type = "salmon",
        tx2gene = transcript_gene_mapping,
        ignoreTxVersion = TRUE,
        ignoreAfterBar = TRUE,
        dropInfReps = FALSE,
        countsFromAbundance = "lengthScaledTPM",
        infRepStat = rowMedians
    )
    return(txi)
}

# Run the benchmark
benchmark_results <- bench::mark(
    tximport_benchmark(),
    filter_gc = FALSE,
    iterations = 11
)

benchmark_results_df <- as.data.frame(
    lapply(benchmark_results, as.character),
    stringsAsFactors = FALSE
)

benchmark_results_df <- benchmark_results_df[, !grepl(
    "result",
    names(benchmark_results_df)
)]

benchmark_results_df <- benchmark_results_df[, !grepl(
    "memory",
    names(benchmark_results_df)
)]

print(benchmark_results)

write.csv(as.data.frame(benchmark_results_df), "tximport_time_memory.csv")
