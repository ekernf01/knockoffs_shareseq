# create local user library path (not present by default)
dir.create(path = Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE)

# Install packages as they were in late 2021
install.packages("versions", lib = Sys.getenv("R_LIBS_USER"))

library(versions)
versions::install.dates(
  pkgs = c("tidyverse", "magrittr", "Matrix", "ggplot2", "irlba", "data.table", "devtools", "Seurat", "BiocManager", "knockoff", "optparse"),
  dates = "2021-11-05",
  lib = Sys.getenv("R_LIBS_USER")
)

BiocManager::install(
  version = 3.14,
  pkgs = c( "DelayedArray", "scater", "scran", "mbkmeans", "HDF5Array" ),
  lib = Sys.getenv("R_LIBS_USER")
)

devtools::install_github("celiaescribe/BDcocolasso", ref = "e883bb34008517adc1bbfb86fbee329b5dc35ab8")
install.packages("~/rlookc", repos = NULL, type = "source", lib = Sys.getenv("R_LIBS_USER"))
