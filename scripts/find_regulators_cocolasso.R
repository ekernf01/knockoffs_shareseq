# Deploy CoCoLASSO on the TF-target inference problem using skin SHARE-seq data.
source("../scripts/setup.R")

# The measurement error model we use is a Gaussian mixture model
# where cells are biologically identical within each cluster
# and all within-cluster variation is due to measurement error.
# Not the only option but a good starting point.
conditions = read.table(
  header = T, 
  text = 
    "knockoff_type condition_on cell_count_cutoff do_downsample
gaussian none  10 FALSE"
)
cat("\nReading in pseudo-bulk data...\n")
normalized_data = set_up_share_skin_pseudobulk(conditions, 1)


# Get RNA data
skin_rna_sce = ???
skin_rna_sce = skin_rna_sce[skin_rna_sce$cluster %in% normalized_data$]
# TODO: Get RNA PC's
# TODO: Get ATAC PC's


cocolasso_vs_chip =function(condition_idx,  conditions){
  attach(conditions[condition_idx,])
  cat("\nRunning condition:\n")
  write.table(conditions[condition_idx,], quote = F, sep = "\t")
  set.seed(seed)
  # Log results
  new_wd = file.path(
    paste0("rna_pca=", rna_pca), 
    paste0("atac_pca=", atac_pca)
  )
  dir.create(new_wd, recursive = T, showWarnings = F)
  setwd(new_wd)
  
  variance_true = normalized_data$pseudo_bulk_rna_var[which(normalized_data$gene_metadata$is_tf),]
  mean_true     = t(normalized_data$pseudo_bulk_rna)[ which(normalized_data$gene_metadata$is_tf),]
  cluster_proportions = normalized_data$pseudo_bulk_metadata$total_cell_count %>% prop.table
  
  DF = 
  # This makes true edges correct; the rest is NA and fails to distinguish between 
  # unknowns and known negatives.
  DF %<>% merge(mouse_chip, all.x = T, all.y = F)
  # This treats every NA as a known negative. 
  DF[["is_verified"]][is.na(DF[["is_verified"]])] = F
  # This fixes unknown edges (those lacking ChIP data). 
  DF[["is_verified"]][ !( DF[["Gene1"]] %in% mouse_chip[["Gene1"]] ) ] = NA
  
  
  DF %>% 
    subset(!is.na(is_verified)) %>%
    subset( q < 0.1, select = c("Gene1", "is_verified")) %>%
    dplyr::group_by(Gene1) %>%
    dplyr::summarise(n_false = sum(!is_verified), n_true = sum(is_verified)) %>%
    merge(gene_metadata) %>% 
    View
  
  # Compute calibration! 
  calibration = data.frame(nominal_fdr = c(1:100)/100, empirical_fdr = NA, num_discoveries = NA)
  i = 0
  for(fdr in calibration$nominal_fdr){
    i = i + 1
    calibration$empirical_fdr[i] =
      DF %>% 
      subset(q<fdr) %>%
      extract2("is_verified") %>% 
      mean(na.rm = T) %>%
      subtract(1, .)
    calibration$num_discoveries[i] = 
      DF %>% 
      subset(q<fdr) %>%
      extract2("is_verified") %>% 
      is.na %>%
      not %>%
      sum
  }
  
  # Add simple binomial standard errors
  calibration %<>% dplyr::mutate( moe_95 = 1.96 * sqrt( empirical_fdr * ( 1 - empirical_fdr ) / num_discoveries ) )
  write.csv(calibration, "calibration/chip_calibration.csv")
  
  # Plot and save calibration
  ggplot(calibration) + 
    geom_point(aes(x = nominal_fdr, y = empirical_fdr)) + 
    geom_errorbar(aes(x = nominal_fdr, 
                      ymin = pmax(0, empirical_fdr - moe_95), 
                      ymax = pmin(1, empirical_fdr + moe_95)
    )) + 
    ggtitle("Calibration vs ChIP-Atlas") + 
    scale_y_continuous(limits = 0:1)
  ggsave("calibration/chip_calibration.pdf", width = 4, height = 4)
  
  setwd(old_wd)
}


conditions = read.table(
  header = T, 
  text = 
    "rna_pca atac_pca seed
       0        0    1
      10        0    1
      20        0    1
      30        0    1
       0       10    1
       0       20    1
       0       30    1
     "
)
parallel::mclapply( seq(nrow(conditions)), cocolasso_vs_chip )

