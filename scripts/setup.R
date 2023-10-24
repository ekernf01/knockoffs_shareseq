# setwd("~/Desktop/jhu/research/projects/knockoffs/applications/share-seq/v13")

# Setup
suppressPackageStartupMessages({
  library("ggplot2")
  library("Matrix")
  library("magrittr")
  library("DelayedArray")
  library("HDF5Array")
  library("scater")
  library("Seurat")
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
create_pseudobulk_path = function(celltype){
  if(  celltype == "keratinocyte" ){
    pseudobulk_path = "input_pseudobulk_keratinocyte"
  } else if(  celltype == "skin" ){
    pseudobulk_path = "input_pseudobulk"
  } else if(  celltype == "pbmc" ){
    pseudobulk_path = "input_pseudobulk_pbmc"
  } else {
    stop(c("celltype must be pbmc, skin, or keratinocyte. Provided value: ", celltype))
  }
  dir.create(pseudobulk_path, showWarnings = F)
  pseudobulk_path
}
# Set up the final results
create_experiment_path = function(conditions, condition_idx){
  attach(conditions[condition_idx,], warn.conflicts = F)
  fp="."
  for(col in colnames(conditions)){
    fp = file.path(
      fp, 
      paste0(
        col, 
        "=",
        conditions[condition_idx,col])
    ) 
  }
  return(fp)
}
atac_peaks = list()

# Load pbmc data
{
  withr::with_dir( file.path(DATALAKE, "multiome_10x", "pbmc"), {
    pbmc_all = HDF5Array::TENxMatrix("10k_PBMC_Multiome_nextgem_Chromium_Controller_filtered_feature_bc_matrix.h5", group = "matrix")
    x = rhdf5::h5read("~/datalake/multiome_10x/pbmc/10k_PBMC_Multiome_nextgem_Chromium_Controller_filtered_feature_bc_matrix.h5", name = "/matrix/features/id")
    rownames(pbmc_all) = as.character(x)
    
    pbmc_rna = pbmc_all[ grepl("ENSG", rownames(pbmc_all)), ]
    pbmc_atac = pbmc_all[!grepl("ENSG", rownames(pbmc_all)), ]
    rownames(pbmc_rna) = as.character(rownames(pbmc_rna))
    rownames(pbmc_atac) = as.character(rownames(pbmc_atac))
    
    # Gene name conversion
    ens_vs_hgnc = biomaRt::getBM(attributes = c("hgnc_symbol", "ensembl_gene_id",  "chromosome_name", "start_position", "end_position"),
                                 filters = "ensembl_gene_id",
                                 values = rownames(pbmc_rna),
                                 mart = biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl"))
    ens_vs_hgnc = setNames(ens_vs_hgnc$hgnc_symbol, nm = ens_vs_hgnc$ensembl_gene_id)
    rownames(pbmc_rna) = ens_vs_hgnc[rownames(pbmc_rna)]
    pbmc_rna = pbmc_rna[!is.na(rownames(pbmc_rna)), ]
    pbmc_rna = pbmc_rna[""!=rownames(pbmc_rna), ]
    pbmc_rna = pbmc_rna[!duplicated(rownames(pbmc_rna)), ]
    
    # Peak filtering
    atac_peaks$pbmc = read.csv("10k_PBMC_Multiome_nextgem_Chromium_Controller_atac_peaks.bed", sep = "", header = F, comment.char = "#")
    peaks_keep = grepl("^chr", atac_peaks$pbmc$V1)
    pbmc_atac = pbmc_atac[peaks_keep,]
    atac_peaks$pbmc = atac_peaks$pbmc[peaks_keep,]
    
    # These cause problems otherwise
    rownames(pbmc_rna) = as.character(rownames(pbmc_rna))
    rownames(pbmc_atac) = as.character(rownames(pbmc_atac))
    
  })
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
    atac_peaks$skin = atac_peaks$keratinocyte = read.csv("GSM4156597_skin.late.anagen.peaks.bed.gz", sep = "\t", header = F)
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
  motif_info$skin <- motifmatchr::matchMotifs(motif_info$all_motifs, GRanges(
    seqnames = atac_peaks$skin$V1, 
    ranges = IRanges(atac_peaks$skin$V2, atac_peaks$skin$V3)), 
    genome = "mm10"
  )
  motif_info$keratinocyte = motif_info$skin
  motif_info$pbmc <- motifmatchr::matchMotifs(
    motif_info$all_motifs,
    GRanges(
      seqnames = atac_peaks$pbmc$V1, 
      ranges = IRanges(atac_peaks$pbmc$V2, atac_peaks$pbmc$V3)
    ), 
    genome = "hg38",
    out = "matches"
  )
}


# Load chip-atlas target gene lists
load_chip_data = function(celltype){
  if(celltype=="keratinocyte"){celltype = "skin"}
  if(celltype=="pbmc"){celltype = "blood"}
  withr::with_dir(
    DATALAKE,
    {
      average_signals = function(DF){
        cat(".")
        new_df = DF[1:2]
        new_df$average_signal = DF[-(1:2)] %>% rowMeans()
        return(new_df)
      }
      chip_files = lapply(
        c("hg19", "hg38", "mm10", "mm9"),
        function(genome) list.files(paste0("chip-atlas/filtered_by_celltype/", celltype, "/", genome), full = T)
      )
      chip_files %<>% Reduce(c, .)
      chip_links =
        lapply(chip_files, function(...) {cat("."); read.csv(...)}, sep = " ", header = T) %>%
        lapply(average_signals) %>%
        lapply(subset, average_signal>mean(average_signal)) %>% # filter for strong signals
        data.table::rbindlist() %>%
        dplyr::mutate(is_verified = T)
    }
  )
  chip_links$Gene1 = chip_links$regulator
  chip_links$Gene2 = chip_links$Target_genes
  chip_links$Gene1 %<>% toupper
  chip_links$Gene2 %<>% toupper
  chip_links$regulator = NULL
  chip_links$Target_genes = NULL
  return (chip_links)
}


get_gene_coords = function(species, gene_symbols){
  gene_coords = NULL
  for(i in 1:100){
    if(!is.null(gene_coords)){break}
    Sys.sleep(1)
    cat(".")
    if(species=="mus_musculus"){
      gene_coords = tryCatch(
        {
          biomaRt::getBM(attributes = c("mgi_symbol", "chromosome_name", "start_position", "end_position"),
                         filters = "mgi_symbol",
                         values = gene_symbols,
                         mart = biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl"))
        }, 
        error = function(e) NULL
      )
    } else if (species=="homo_sapiens") {
      gene_coords = tryCatch(
        {
          biomaRt::getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
                         filters = "hgnc_symbol",
                         values = gene_symbols,
                         mart = biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl"))
        }, 
        error = function(e) NULL
      )
    } else {
      stop(paste0("species must be 'homo_sapiens' or 'mus_musculus'. Provided: ",  species))
    }
  }
  return(gene_coords)
}

#' Obtain a subset of TF-target hypotheses supported by a motif in a region whose chromatin accessibility
#' correlates with the target gene expression.
#' 
get_motif_supported_hypotheses = function(normalized_data, celltype){
  cat("Getting gene coordinates from ensembl. Each dot is one attempt (the download sometimes fails).\n")
  gene_coords = get_gene_coords(species = normalized_data$species, gene_symbols = normalized_data$gene_metadata$Gene1)
  gene_coords = gene_coords[!duplicated(gene_coords[[1]]),]
  stopifnot(length(gene_coords[[1]])==length(normalized_data$gene_metadata$Gene1))
  gene_coords = GRanges(
    seqnames = paste0("chr", gene_coords$chromosome_name),
    ranges = IRanges(gene_coords$start_position, gene_coords$end_position),
    genome = ifelse(normalized_data$species=="mus_musculus", "mm10", "hg38"), 
    symbol = gene_coords[[1]]
  )
  atac_peaks_gr = GRanges(
    seqnames = atac_peaks[[celltype]]$V1,
    ranges = IRanges(atac_peaks[[celltype]]$V2, atac_peaks[[celltype]]$V3),
    genome = ifelse(normalized_data$species=="mus_musculus", "mm10", "hg38")
  )
  overlaps = GenomicRanges::findOverlaps(
    gene_coords,
    atac_peaks_gr,
    maxgap = 1e5
  )
  overlaps %<>% as.data.frame() %>% set_colnames(c("gene", "enhancer"))
  overlaps[["distance"]] = GenomicRanges::distance(
    gene_coords[overlaps[["gene"]]],
    atac_peaks_gr[overlaps[["enhancer"]]]
  )
  overlaps$correlation = NA
  cat("Correlating genes with nearby ATAC peaks.\n")
  for(i in seq_along(overlaps[[1]])){
    if((i%%10000)==0){cat(i, " of ", length(overlaps[[1]]), "\n")}
    peak_idx = overlaps[i, "enhancer"]
    gene_idx = gene_coords[overlaps[i, "gene"]]$symbol
    overlaps[i, "correlation"] = cor(
      normalized_data$pseudo_bulk_atac[,peak_idx],
      normalized_data$pseudo_bulk_rna[,gene_idx],
    )
  }
  overlaps %<>% dplyr::mutate(is_kept = correlation > 0 | distance < 2e3)
  ggplot2::ggplot(overlaps) +
    geom_hex(aes(distance, correlation)) +
    facet_wrap(~is_kept)
  ggtitle("Enhancer-gene pairing in skin SHARE-seq data")
  ggsave("enhancer_pairing.pdf", width = 5, height = 5)
  enhancer_gene_links = subset(overlaps, is_kept)
  enhancer_gene_links$Gene2 = colnames(normalized_data$pseudo_bulk_rna)[enhancer_gene_links$gene]
  motif_enhancer_links = 
    motif_info[[celltype]] %>% 
    assay("motifMatches") %>% 
    summary %>% 
    set_colnames(c("enhancer", "motif_index", "is_connected"))
  gene_motif_links = link_genes_to_motifs(motif_info)
  motif_gene_links = merge(motif_enhancer_links, enhancer_gene_links, by = "enhancer")
  gene_gene_links_from_motif_analysis = merge(gene_motif_links, motif_gene_links, by = "motif_index")
  gene_gene_links_from_motif_analysis = gene_gene_links_from_motif_analysis[c("Gene1", "Gene2")]
  return(gene_gene_links_from_motif_analysis)
}

link_genes_to_motifs = function(motif_info){
  motif_info$all_motifs %>%
    sapply(TFBSTools::name) %>%
    data.frame(genes=., motif = names(.)) %>%
    (dplyr::add_rownames) %>%
    dplyr::rename(motif_index = rowname) %>%
    dplyr::mutate(genes = gsub("\\(var\\..*\\)", "", genes)) %>%
    tidyr::separate_rows("genes", sep = "::") %>%
    dplyr::mutate(Gene1 = toupper(genes)) %>%
    extract(c("motif_index", "motif", "Gene1"))
}

#' Load the full data in a nice format
#' 
#' list with one SingleCellExperiment each for RNA and ATAC.
#' RNA are normalized by total counts.
#' the RNA cell barcodes should be equal to, or a subset of, the ATAC ones. 
#' 
set_up_sce = function(celltype = "skin"){
  if(celltype=="pbmc"){
    rownames(pbmc_rna) = as.character(rownames(pbmc_rna))
    rownames(pbmc_atac) = as.character(rownames(pbmc_atac))
    rna_sce  = SingleCellExperiment(assays = list(counts = pbmc_rna))
    atac_sce = SingleCellExperiment(assays = list(counts = pbmc_atac))
    colData( rna_sce)["atac.bc"] = colnames(pbmc_rna)
    colData(atac_sce)["atac.bc"] = colnames(pbmc_rna)
    colData(rna_sce)[["celltype"]] = "unknown"
  } else if(celltype %in% c("skin", "keratinocyte")){
    rna_sce = SingleCellExperiment(assays = list(counts = t(skin_rna)))
    colData(rna_sce)[["rna.bc"]] = rownames(skin_rna)
    colData(rna_sce) = merge(colData(rna_sce),
                             share_seq_skin_metadata,
                             by = "rna.bc",
                             all.x = T,
                             all.y = F)
    # Exclude cells with no ATAC counterpart or cell type label.
    rna_sce = rna_sce[,!is.na(rna_sce[["celltype"]])]
    
    # Set up paired ATAC data
    atac_sce = SingleCellExperiment(assays = list(counts = skin_atac))
    atac_sce
    colData(atac_sce)[["atac.bc"]] = rownames(colData(atac_sce))
    colData(atac_sce) = merge(colData(atac_sce),
                              colData(rna_sce),
                              by = "atac.bc",
                              all.x = T,
                              all.y = F)
  } else {
    stop(c("celltype must be pbmc, skin, or keratinocyte. Provided value: ", celltype))
  }
  # Exclude cells we don't have ChIP data for
  if(celltype == "keratinocyte"){
    keep = c("keratinocyte",
             "hair_shaft",
             "bulge")
    rna_sce = rna_sce[,rna_sce[["supertype"]] %in% keep]
  }
  
  # Normalize
  sizeFactors(rna_sce) <- colSums(assay(rna_sce, "counts"))
  rna_sce <- scater::logNormCounts(rna_sce)
  
  return(list( rna_sce = rna_sce,
               atac_sce = atac_sce ))
}

# Reload cluster average expression.
# Only works after cluster_cells.R has been run.
set_up_share_skin_pseudobulk = function(conditions, i){
  # Which experiment are we doing?
  celltype         = conditions$celltype[[i]]
  knockoff_type     = conditions$knockoff_type[[i]]
  error_mode        = conditions$error_mode[[i]]
  cell_count_cutoff = conditions$cell_count_cutoff[[i]]
  condition_on      = conditions$condition_on[[i]]
  
  # Get TF's
  if(celltype %in% c("skin", "keratinocyte")){
    species = "mus_musculus"
    withr::with_dir(
      DATALAKE,
      {  
        tf_fields = c("Species",   "Symbol" ,   "Ensembl", "Family", "Entrez.ID")
        tf_list = rbind(
          read.csv("mouse_tfs/Mus_musculus_TF.txt", sep = "\t")[tf_fields],
          read.csv("mouse_tfs/Mus_musculus_TF_cofactors.txt", sep = "\t")[tf_fields]
        )
      })
  } else if(celltype =="pbmc"){
    species = "homo_sapiens"
    withr::with_dir(
      DATALAKE,
      { 
        tf_list = read.csv("human_tfs/humanTFs.csv", row.names = 1)[1:6] %>% subset(Is.TF.=="Yes")
        colnames(tf_list)[1:2] = c( "Ensembl" , "Symbol"  )
      }
    )
  } else {
    stop(paste0("celltype must be 'skin', 'keratinocyte', or 'pbmc'; got ", celltype))
  }
  
  
  # Load expression pseudo-bulk data
  withr::with_dir(create_pseudobulk_path(celltype), {
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
    is_tf = rownames(pseudo_bulk_rna) %in% tf_list$Symbol
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
  pseudo_bulk_rna_var %<>% t
  pseudo_bulk_rna_noisy %<>% t
  pseudo_bulk_atac %<>% t
  return(list(
    non_tf_expression    = pseudo_bulk_rna[,!gene_metadata$is_tf],
    tf_expression        = pseudo_bulk_rna[, gene_metadata$is_tf],
    tf_expression_noisy  = pseudo_bulk_rna_noisy[, gene_metadata$is_tf],
    gene_metadata              = gene_metadata,
    pseudo_bulk_metadata       = pseudo_bulk_metadata,
    pseudo_bulk_rna            = pseudo_bulk_rna,
    pseudo_bulk_rna_noisy      = pseudo_bulk_rna_noisy,
    pseudo_bulk_rna_var        = pseudo_bulk_rna_var,
    pseudo_bulk_atac           = pseudo_bulk_atac,
    motif_activity             = ( pseudo_bulk_atac %*% motifmatchr::motifMatches(motif_info[[celltype]]) ) , 
    species = species
  ) )
}

# This is a fast way to compute variable importance statistics.
dfmax = 21
fast_lasso_penalty = function(X, X_k, y) {
  cat(".")
  suppressWarnings(
    knockoff::stat.lasso_lambdasmax(
      y   = y,
      X   = X,
      X_k = X_k,
      dfmax = dfmax
    )
  )
}

#' Z-score the data, adding noise to constants to avoid outputting anything with zero sd
#'
safe_scale = function(x, ...){
  if(sd(x)==0){
    return( x + rnorm(length(x)) )
  } else {
    return( (x - mean(x)) / sd(x) )
  }
}

