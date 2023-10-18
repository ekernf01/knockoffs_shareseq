# create local user library path (not present by default)
dir.create(path = Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE, mode = "777")
.libPaths(Sys.getenv("R_LIBS_USER"))

# # Install packages as they were in late 2021
install.packages("https://cran.r-project.org/src/contrib/Archive/remotes/remotes_2.4.2.tar.gz",
                 source = TRUE,
                 repos = NULL)
pv  = read.table(header = 1, text =
  "package version
tidyverse 1.3.1
magrittr 2.0.3
Matrix 1.6.1
ggplot2 3.4.2
irlba 2.3.3
data.table 1.14.8
devtools 2.4.2
Seurat 4.30.0.1
BiocManager 1.30.16
knockoff 0.3.3
vita 1.0.0
optparse 1.7.1"
)
for(i in rownames(pv)){
  remotes::install_version(
    pv[i, "package"],
    version = pv[i, "version"],
    upgrade = "never",
    quiet = TRUE,
    repos = "https://cloud.r-project.org"
  )
}

BiocManager::install(
  version = "3.14",
  pkgs = c(
    "DelayedArray",
    "scater",
    "scran",
    "mbkmeans",
    "JASPAR2018",
    "HDF5Array",
    "motifmatchr",
    'BiocGenerics',
    'DelayedMatrixStats',
    'biomaRt',
    "GenomicRanges",
    "BSgenome.Hsapiens.UCSC.hg38",
    "BSgenome.Mmusculus.UCSC.mm10"
  )
)


# Install our package
remotes::install_github("ekernf01/rlookc", ref="840bf3a")
remotes::install_github("ekernf01/gm", ref="547a42a")

