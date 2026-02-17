# List of CRAN packages
cran_packages <- c(
  "readr", "dplyr", "data.table", "parallel", "caret", "glmnet",
  "paletteer", "stringr", "ComplexHeatmap", "tidyr", "ggbreak",
  "ggplot2", "ggridges", "tibble", "ggpubr", "circlize", "ggrepel",
  "gridExtra", "RMariaDB", "patchwork", "purrr", "reticulate",
  "scales", "igraph", "ggraph", "tidygraph", "ggbeeswarm",
  "ggseqlogo", "tools", "ggrastr", "Cairo", "fmsb", "ggsci", "colorRamp2"
)

# List of Bioconductor packages
bioc_packages <- c(
  "clusterProfiler", "org.Hs.eg.db", "GenomicRanges", "rtracklayer",
  "IlluminaHumanMethylationEPICv2anno.20a1.hg38",
  "TxDb.Hsapiens.UCSC.hg38.knownGene", "DESeq2",
  "GenomicFeatures", "biomaRt", "methylCIPHER", "UpSetR",
  "ChAMP", "GO.db", "Biostrings", "DunedinPACE"
)

# Install missing CRAN packages
install_if_missing_cran <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org/")
  }
}
invisible(lapply(cran_packages, install_if_missing_cran))

# Install missing Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org/")
}

install_if_missing_bioc <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}
invisible(lapply(bioc_packages, install_if_missing_bioc))

# Load all packages
all_packages <- c(cran_packages, bioc_packages)
invisible(lapply(all_packages, library, character.only = TRUE))

message("All packages are installed and loaded successfully!")
