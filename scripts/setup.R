# setwd("~/Desktop/jhu/research/projects/knockoffs/applications/share-seq/v13")

# Setup
suppressPackageStartupMessages({
  library("ggplot2")
  library("Matrix")
  library("magrittr")
  library("DelayedArray")
  library("HDF5Array")
  library("scater")
  library("scran")
  library("optparse")
})
set.seed(0)

# for eric laptop
DATALAKE = "~/Desktop/jhu/research/datalake"
# for AWS
if(!dir.exists(DATALAKE)){
  DATALAKE = "~/datalake"
}
cat("Looking for data in", DATALAKE, "\n")
# otherwise, complain.
if(!dir.exists(DATALAKE)){
  stop("Datalake not found. Place it in '~/datalake' or Modify `setup.R`.\n")
}

# Set up some key intermediate steps
create_pseudobulk_path = function(keratinocyte_only){
  if(  keratinocyte_only ){
    pseudobulk_path = "input_pseudobulk_keratinocyte"
  } else {
    pseudobulk_path = "input_pseudobulk"
  }
  dir.create(pseudobulk_path, showWarnings = F)
  pseudobulk_path
}
# Set up the final results
create_experiment_path = function(conditions, condition_idx){
  attach(conditions[condition_idx,], warn.conflicts = F)
  file.path(
    paste0("keratinocyte_only=", keratinocyte_only),
    paste0("cell_count_cutoff=", cell_count_cutoff),
    paste0("condition_on=", condition_on),
    paste0("knockoff_type=", knockoff_type),
    paste0("error_mode=", error_mode)
  )
}

# Load skin data
{
  withr::with_dir( file.path(DATALAKE, "share_seq", "skin"), {
    skin_rna  = HDF5Array::TENxMatrix("GSM4156608_skin.late.anagen.rna.counts.h5", "RNA")
    skin_atac = HDF5Array::TENxMatrix("GSM4156597_skin.late.anagen.counts.atac.h5", "ATAC")
    share_seq_skin_metadata = read.csv("GSM4156597_skin_celltype.txt.gz",
                                       sep ="\t")
    colnames(skin_rna) = read.csv("GSM4156608_skin.late.anagen.rna.counts.txt.gz.rownames.txt",
                                  header = T)[[1]]
    rownames(skin_rna) = colnames(read.csv("GSM4156608_skin.late.anagen.rna.counts.txt.gz",
                                           nrow = 1, sep = "\t", row.names = 1,
                                           header = T))
    skin_atac_peaks = read.csv("GSM4156597_skin.late.anagen.peaks.bed.gz", sep = "\t", header = F)
    skin_atac_bc    = read.csv("GSM4156597_skin.late.anagen.barcodes.txt.gz", header = F)[[1]]
  })
  # The last barcode is off by 48.
  length(intersect(skin_atac_bc, share_seq_skin_metadata$atac.bc))
  get_final_barcode = function(x) x %>% strsplit("\\.") %>% lapply(rev) %>% sapply(extract2, 1) %>% as.numeric
  corrected_final_barcode = skin_atac_bc %>% get_final_barcode %>% magrittr::subtract(48) %>% formatC(width = 2, format = "d", flag = "0")
  table(corrected_final_barcode)
  skin_atac_bc = gsub("[0-9]{2}$", "", skin_atac_bc)
  skin_atac_bc = paste0(skin_atac_bc, corrected_final_barcode)
  colnames(skin_atac) = skin_atac_bc
  length(intersect(skin_atac_bc, share_seq_skin_metadata$atac.bc))

  # Include cell type annotations
  
  supertypes = read.csv(text =
                          "celltype,supertype
ahighCD34+ bulge,bulge
alowCD34+ bulge,bulge
Basal,keratinocyte
Dermal Fibroblast,other
Dermal Papilla,other
Dermal Sheath,other
Endothelial,other
Hair Shaft-cuticle.cortex,hair_shaft
Infundibulum,hair_shaft
IRS,hair_shaft
K6+ Bulge Companion Layer,hair_shaft
Medulla,hair_shaft
Melanocyte,melanocyte
Mix,other
ORS,hair_shaft
Sebaceous Gland,other
Spinous,keratinocyte
TAC-1,keratinocyte
TAC-2,keratinocyte"
  )
  share_seq_skin_metadata %<>% merge(supertypes, by = "celltype")
}

# Skin ATAC feature extraction & pruning: global motif activity per cell and motif-peak-gene network construction
{  
  motif_info = list(
    mouse_motifs <- TFBSTools::getMatrixSet(
      JASPAR2018::JASPAR2018, 
      opts = list(
        "species" = c(10090), # mouse
        "all_versions" = TRUE
      )
    ),
    human_motifs <- TFBSTools::getMatrixSet(
      JASPAR2018::JASPAR2018, 
      opts = list(
        "species" = c( 9606), # human
        "all_versions" = TRUE
      )
    ),
    all_motifs = c(mouse_motifs, human_motifs)
  )
  motif_info$motif_ix <- motifmatchr::matchMotifs(motif_info$all_motifs, GRanges(
    seqnames = skin_atac_peaks$V1, 
    ranges = IRanges(skin_atac_peaks$V2, skin_atac_peaks$V3)), 
    genome = "mm10"
  )
}


# Load chip-atlas target gene lists
{
  cat("\nLoading ChIP data...\n")
  average_signals = function(DF){
    new_df = DF[1:2]
    new_df$average_signal = DF[-(1:2)] %>% rowMeans()
    return(new_df)
  }
  withr::with_dir(
    DATALAKE,
    {
      chip_files = list.files(
        paste0("chip-atlas/filtered_by_celltype/skin/", 
               c("hg19",
                 "hg38",
                 "mm10",
                 "mm9")
        ), 
        full = T)
      mouse_chip =
        lapply(chip_files, read.table, header = T) %>%
        lapply(average_signals) %>%
        lapply(subset, average_signal>mean(average_signal)) %>% 
        data.table::rbindlist() %>%
        dplyr::mutate(is_verified = T)
    }
  )
  frequency_in_chip_data =
    mouse_chip[["regulator"]] %>%
    table %>%
    as.data.frame() %>%
    set_colnames(c("Gene1", "chip_freq"))
}

# Load the full data in a nice format
set_up_sce = function(keratinocyte_only = T){
  skin_rna_sce = SingleCellExperiment(assays = list(counts = t(skin_rna)))
  colData(skin_rna_sce)[["rna.bc"]] = rownames(skin_rna)
  colData(skin_rna_sce) = merge(colData(skin_rna_sce),
                                share_seq_skin_metadata,
                                by = "rna.bc",
                                all.x = T,
                                all.y = F)
  # Exclude cells with no ATAC counterpart or cell type label.
  skin_rna_sce = skin_rna_sce[,!is.na(skin_rna_sce[["celltype"]])]

  # Exclude cells we don't have ChIP data for
  if(keratinocyte_only){
    keep = c("keratinocyte",
      "hair_shaft",
      "bulge")
    skin_rna_sce = skin_rna_sce[,skin_rna_sce[["supertype"]] %in% keep]
  }

  # Normalize
  sizeFactors(skin_rna_sce) <- colSums(assay(skin_rna_sce, "counts"))
  skin_rna_sce <- scater::logNormCounts(skin_rna_sce)

  # Retrieve metadata, tSNE, etc, if it exists
  withr::with_dir(create_pseudobulk_path(keratinocyte_only), {
    cell_metadata = NULL
    try({cell_metadata = read.csv("rna_metadata.csv")}, silent = T)
  })

  # Set up paired ATAC data
  skin_atac_sce = SingleCellExperiment(list("counts" = skin_atac))
  skin_atac_sce
  colData(skin_atac_sce)[["atac.bc"]] = rownames(colData(skin_atac_sce))
  colData(skin_atac_sce) = merge(colData(skin_atac_sce),
                                 colData(skin_rna_sce),
                                 by = "atac.bc",
                                 all.x = T,
                                 all.y = F)
  return(list( skin_rna_sce = skin_rna_sce,
               skin_atac_sce = skin_atac_sce,
               cell_metadata = cell_metadata ))
}

# Reload cluster average expression.
# Only works after cluster_cells.R has been run.
set_up_share_skin_pseudobulk = function(conditions, i){
  # Which experiment are we doing?
  keratinocyte_only = conditions$keratinocyte_only[[i]]
  knockoff_type     = conditions$knockoff_type[[i]]
  error_mode        = conditions$error_mode[[i]]
  cell_count_cutoff = conditions$cell_count_cutoff[[i]]
  condition_on      = conditions$condition_on[[i]]

  # Get TF's
  tf_fields = c("Species",   "Symbol" ,   "Ensembl", "Family", "Entrez.ID")
  withr::with_dir(
    DATALAKE,
    {
      mouse_tf_atfdb = rbind(
        read.csv("mouse_tfs/Mus_musculus_TF.txt", sep = "\t")[tf_fields],
        read.csv("mouse_tfs/Mus_musculus_TF_cofactors.txt", sep = "\t")[tf_fields]
      )
    }
  )

  # Load expression pseudo-bulk data
  withr::with_dir(create_pseudobulk_path(keratinocyte_only), {
    pseudo_bulk_atac = NULL
    try({pseudo_bulk_atac = read.csv("atac_pseudo_bulk.csv",     row.names = 1) %>% as.matrix})
    pseudo_bulk_rna       = read.csv("rna_pseudo_bulk.csv",      row.names = 1) %>% as.matrix
    pseudo_bulk_rna_var   = read.csv("rna_pseudo_bulk_vars.csv", row.names = 1) %>% as.matrix
    pseudo_bulk_metadata  = read.csv("metacell_metadata.csv",    row.names = 1)
    non_empty_clusters = colnames(pseudo_bulk_rna) %>% gsub("^X", "", .) %>% as.numeric
    pseudo_bulk_metadata %<>% subset(cluster %in% non_empty_clusters)
    pseudo_bulk_rna_var = pseudo_bulk_rna_var[,non_empty_clusters]
  })

  # Exclude clusters with too few cells.
  write.table( table(pseudo_bulk_metadata$total_cell_count >= cell_count_cutoff), "n_clusters_passing_cutoff.csv" )
  metacells_keep = pseudo_bulk_metadata$total_cell_count >= cell_count_cutoff
  try({pseudo_bulk_atac     %<>% extract( , metacells_keep)})
  pseudo_bulk_rna      %<>% extract( , metacells_keep)
  pseudo_bulk_rna_var  %<>% extract( , metacells_keep)
  pseudo_bulk_metadata %<>% extract( metacells_keep , )
  withr::with_dir( DATALAKE, {
    skin_atac_peaks = read.csv("share_seq/skin/GSM4156597_skin.late.anagen.peaks.bed.gz",
                               sep = "\t", header = F)
  })
  dim(mouse_tf_atfdb)

  # Downsample or resample to study effect of measurement error
  pseudo_bulk_rna_noisy = pseudo_bulk_rna
  if( tolower(error_mode) == "downsample" ){
    pseudo_bulk_rna_noisy = Seurat::SampleUMI(pseudo_bulk_rna, max.umi = 1e4) %>% as.matrix
  }
  if( tolower(error_mode) == "resample" ){
    for(  i in seq(nrow(pseudo_bulk_rna))){
      for(j in seq(ncol(pseudo_bulk_rna))){
        pseudo_bulk_rna_noisy[i,j] = rpois(1, pseudo_bulk_rna[i, j])
      }
    }
  }

  # Normalize by total counts. Compute some gene-level summaries.
  pseudo_bulk_atac      = sweep(pseudo_bulk_atac,     2, colSums(pseudo_bulk_atac      ), FUN = "/")*1e6
  pseudo_bulk_rna       = sweep(pseudo_bulk_rna,       2, colSums(pseudo_bulk_rna      ), FUN = "/")*1e6
  pseudo_bulk_rna_noisy = sweep(pseudo_bulk_rna_noisy, 2, colSums(pseudo_bulk_rna_noisy), FUN = "/")*1e6
  gene_metadata = data.frame(
    Gene1 = rownames(pseudo_bulk_rna),
    mean_expression = rowMeans(pseudo_bulk_rna),
    sd_expression   = apply(   pseudo_bulk_rna, 1, sd),
    is_tf = rownames(pseudo_bulk_rna) %in% mouse_tf_atfdb$Symbol
  )
  gene_metadata$cv_expression = with(gene_metadata, sd_expression / mean_expression)

  # Exclude genes below 1 CPM (mean).
  genes_keep = gene_metadata$mean_expression >= 1
  pseudo_bulk_rna       = pseudo_bulk_rna[       genes_keep, ]
  pseudo_bulk_rna_noisy = pseudo_bulk_rna_noisy[ genes_keep, ]
  pseudo_bulk_rna_var   = pseudo_bulk_rna_var[   genes_keep, ]
  gene_metadata         =   gene_metadata[genes_keep, ]

  # Knockoff filter code expects 1 variable per column, not per row
  pseudo_bulk_rna %<>% t
  pseudo_bulk_rna_noisy %<>% t
  return(list(
    mouse_non_tf_expression    = pseudo_bulk_rna[,!gene_metadata$is_tf],
    mouse_tf_expression        = pseudo_bulk_rna[, gene_metadata$is_tf],
    mouse_tf_expression_noisy  = pseudo_bulk_rna_noisy[, gene_metadata$is_tf],
    gene_metadata              = gene_metadata,
    pseudo_bulk_metadata       = pseudo_bulk_metadata,
    pseudo_bulk_rna            = pseudo_bulk_rna,
    pseudo_bulk_rna_noisy      = pseudo_bulk_rna_noisy,
    pseudo_bulk_rna_var        = pseudo_bulk_rna_var,
    pseudo_bulk_atac           = pseudo_bulk_atac,
    motif_activity = t(motifmatchr::motifMatches(motif_info$motif_ix)) %*% pseudo_bulk_atac
  ) )
}
