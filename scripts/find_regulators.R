# Deploy knockoff-based variable selection
# on the TF-target inference problem.
source("../scripts/setup.R")
conditions = read.table(
  header = T,
  text =
    "knockoff_type tf_activity_type condition_on cell_count_cutoff error_mode seed cell_types require_motif_support only_motif_support
gaussian   rna none  10       none 1 skin F  F
gaussian   rna none 100       none 1 skin F  F
gaussian   rna none 500       none 1 skin F  F
permuted   rna none  10       none 1 skin F  F
permuted   rna none 100       none 1 skin F  F
permuted   rna none 500       none 1 skin F  F
gaussian   rna none 500 downsample 1 skin F  F
gaussian   rna none  10   resample 1 skin F  F
gaussian   rna none 100   resample 1 skin F  F
gaussian   rna none 500   resample 1 skin F  F
gaussian motif none  10       none 1 skin F  F
gaussian motif none 100       none 1 skin F  F
gaussian motif none 500       none 1 skin F  F
gaussian  both none  10       none 1 skin F  F
gaussian  both none 100       none 1 skin F  F
gaussian  both none 500       none 1 skin F  F
gaussian  both  pca  10       none 1 skin F  F
gaussian  both  pca 100       none 1 skin F  F
gaussian  both  pca 500       none 1 skin F  F
gaussian  both  pca  10       none 1 skin  T F
gaussian  both  pca 100       none 1 skin  T F
gaussian  both  pca 500       none 1 skin  T F
gaussian  both  pca  10       none 1 skin F   T
gaussian  both  pca 100       none 1 skin F   T
gaussian  both  pca 500       none 1 skin F   T
")
conditions = Reduce(
  rbind, 
  list(
    conditions, 
    dplyr::mutate(conditions, cell_types="keratinocyte"),
    dplyr::mutate(conditions, cell_types="pbmc")
  )
)

write.csv(conditions, "experiments_to_run.csv")
old_wd = getwd()
condition_idx = 3 # For interactive debugging

#' Run a single experiment described by a row of the 'conditions' dataframe defined above.
#'
#' @param knockoff_type 'permuted' or 'gaussian'. 
#' @param tf_activity 'rna' or 'motif' or 'both'. How to extract TF activities. 
#' @param condition_on 'none' or 'pca.' if 'pca', then include the top 10 PC's from both ATAC and RNA data.
#' @param cell_count_cutoff If a cluster has fewer cells than this, it gets excluded from our analysis.
#' @param error_mode If 'none', take no special action. If 'resample', replace mRNA count Xij with a Poisson(Xij). If 'downsample', use Seurat's SampleUMI function to thin out the data.
#' @param seed Random seed in case you want to replicate experiments. 
#' @param cell_type skin: all the mouse skin SHARE data. keratinocyte: a subset of the mouse skin SHARE data. pbmc: the human data. 
#' @param require_motif_support If T, retain an edge only when there is a motif for the TF in a nearby open chromatin region.
#' @param only_motif_support If T, ignore everything else and just get the FDR you'd see if you used motif-matching as the sole inference method.
#' 
do_one = function(condition_idx, reuse_results = F, spread_load = T){
  set.seed(conditions[condition_idx,"seed"])
  if(spread_load){  
    Sys.sleep(120*(condition_idx-1)) #spread out peak RAM
  }
  dir.create("logs")
  withr::with_output_sink(file.path("logs", condition_idx), {
    withr::with_message_sink(file.path("logs", condition_idx), {
      # Set up environment
      attach(conditions[condition_idx,], warn.conflicts = F)
      cat("\n", as.character(Sys.time()), "\n")
      cat("\nRunning condition:\n")
      write.table(t(conditions[condition_idx,]), quote = F, sep = "\t")
      if(!reuse_results){
        normalized_data = set_up_share_skin_pseudobulk(conditions, condition_idx)
      }
      new_wd = create_experiment_path(conditions, condition_idx)
      dir.create(new_wd, recursive = T, showWarnings = F)
      withr::with_dir( new_wd, {
        if(!only_motif_support & !reuse_results){

          # For this analysis, we want expression centered and scaled
          # For genes with zero variance, esp after downsampling,
          # there's too much brittle software downstream. (Not all of it mine!)
          # So we randomize instead of keeping constants or zeroes.
          cat("\n", as.character(Sys.time()), "\n")
          cat("\nScaling RNA data...\n")
          
          normalized_data$mouse_non_tf_expression   %<>% apply(2, safe_scale)
          normalized_data$mouse_tf_expression       %<>% apply(2, safe_scale)
          normalized_data$mouse_tf_expression_noisy %<>% apply(2, safe_scale)
          normalized_data$motif_activity %<>% apply(2, safe_scale)
          # Choose a measure of TF activity
          if(tf_activity_type == "motif"){
            tf_activity = normalized_data$motif_activity
            colnames(tf_activity) = colnames(normalized_data$motif_activity)
          } else if(tf_activity_type=="both"){
            tf_activity =       cbind(         normalized_data$motif_activity,           normalized_data$mouse_tf_expression)
            colnames(tf_activity) = c(colnames(normalized_data$motif_activity), colnames(normalized_data$mouse_tf_expression))
          } else if(tf_activity_type=="rna"){
            tf_activity = normalized_data$mouse_tf_expression
            colnames(tf_activity) = colnames(normalized_data$mouse_tf_expression)
          } else{
            stop("tf_activity_type must be 'rna' or 'motif' or 'both'")
          }
          # Add simulated measurement error to TF activity
          if(error_mode == "none"){
            tf_activity_before_noise = tf_activity
          } else {
            stopifnot("Unless tf_activity_type is 'rna', error_mode must be 'none'."=tf_activity_type=="rna")
            tf_activity = normalized_data$mouse_tf_expression_noisy
            tf_activity_before_noise = normalized_data$mouse_tf_expression
          } 
          # Add covariates / surrogate variables
          if(condition_on == "pca"){
            atac_pca = irlba::prcomp_irlba(normalized_data$pseudo_bulk_atac, n = 10)
            rna_pca = irlba::prcomp_irlba(normalized_data$pseudo_bulk_rna, n = 10)
            covariates = cbind( tf_activity, atac_pca$x, rna_pca$x)
          } else if(condition_on == "none"){
            covariates = tf_activity
          } else {
            stop("condition_on must be 'pca' or 'none'.")
          }
          # Generate knockoffs
          cat("\n", as.character(Sys.time()), "\n")
          cat("\nGenerating knockoffs\n")
          if(knockoff_type == "permuted"){
            # Assumes variables are in columns! Otherwise you're
            # scrambling all genes within a sample, not all samples within a gene.
            knockoffs = apply(covariates, 2, sample, replace = F)
          } else if (knockoff_type=="gaussian"){
            shrink_by = corpcor::estimate.lambda(covariates)
            Sigma = corpcor::cov.shrink(x = covariates, lambda = shrink_by, lambda.var = 0 )
            Sigma = as(Sigma, "matrix")
            knockoffs = rlookc::createHighDimensionalKnockoffs(
              covariates, lambda = shrink_by,
            )
          } else {
            stop("Unknown type of knockoffs requested by upstream code\n.")
          }
          # Strip out extra covariates
          knockoffs = knockoffs[,1:ncol(tf_activity)]
          # Save knockoff realization
          dir.create("output_knockoffs", recursive = T, showWarnings = F)
          saveRDS(knockoffs, "output_knockoffs/knockoffs.Rda")
          saveRDS(tf_activity, "output_knockoffs/original_features.Rda")
          
          # Check calibration with simulated targets
          cat("\n", as.character(Sys.time()), "\n")
          cat("\nChecking calibration with simulated targets\n")
          dir.create("calibration", showWarnings = F, recursive = T)
          calibration_typical = rlookc::calibrate__simulateY(
            n_sim = 100,
            active_set_size = pmax(1, rpois(2, n = 100)),
            X = tf_activity_before_noise,
            X_observed = tf_activity,
            knockoffs = knockoffs,
            statistic = fast_lasso_penalty, 
            plot_savepath = "calibration/average_case_calibration.pdf"
          )
          write.csv(calibration_typical$calibration$fdr, "calibration/average_case_calibration.csv")
          
          
          # Find regulators of only non-TF genes.
          # We do this to avoid "spousal problems", since spouses may be linked in
          # the MRF structure needed to express a given causal DAG
          # even when they are not linked in the corresponding DAG.
          cat("\n", as.character(Sys.time()), "\n")
          cat("\nChecking calibration with real targets\n")
          cat("Computing knockoff statistics. This may take a while.\n")
          w = lapply(
            seq(ncol(normalized_data$mouse_non_tf_expression)),
            function(i){
              if(i %% 100 == 0){ cat( "\n", i, "/", ncol(normalized_data$mouse_non_tf_expression) )}
              fast_lasso_penalty(
                y   = normalized_data$mouse_non_tf_expression[,i],
                X   = tf_activity,
                X_k = knockoffs
              )
            }
          )
          length(w[[1]])
          length(w)
          cat("Saving knockoff statistics.\n")
          saveRDS(w, "output_knockoffs/w.Rds")
          
          # Assemble results nicely
          cat("\n", as.character(Sys.time()), "\n")
          cat("Assembling knockoff stats and cleaning up.\n")
          DF_tf_target = list()
          for(k in seq(ncol(normalized_data$mouse_non_tf_expression))){
            if(k %% 100 == 0){ cat( "\n", k, "/", ncol(normalized_data$mouse_non_tf_expression) )}
            DF_tf_target[[k]] = data.frame(
              Gene1 = colnames(tf_activity),
              Gene2 = colnames(normalized_data$mouse_non_tf_expression)[ k],
              knockoff_stat = w[[k]]
            ) %>%
              subset(abs(knockoff_stat) > 0)
          }
          rm("w"); gc()
          DF_tf_target = data.table::rbindlist(DF_tf_target)
          
          # Save results
          write.csv(DF_tf_target, "output_knockoffs/knockoff_stats.csv")
          rm("DF_tf_target"); gc()
        }
        
        # We will try using 1) all testable hypotheses, with 2) testable hypotheses motif support, or 3) as a control, only motif-based hypotheses
        cat("\n", as.character(Sys.time()), "\n")
        cat("\nIntegrating knockoff-based and motif-based hypotheses... \n")
        if( only_motif_support ){
          DF = get_motif_supported_hypotheses(normalized_data)
          DF$Gene1 %<>% toupper()
          DF$Gene2 %<>% toupper()
          DF$knockoff_stat = 100 + rnorm(nrow(DF)) #these fake knockoff statistics will yield q=0 for all hypotheses
        } else {
          DF = read.csv("output_knockoffs/knockoff_stats.csv")
          DF$Gene1 %<>% toupper()
          DF$Gene2 %<>% toupper()
          if( require_motif_support ){
            mosh = get_motif_supported_hypotheses(normalized_data)
            mosh$Gene1 %<>% toupper()
            mosh$Gene2 %<>% toupper()
            DF = merge(DF, mosh, type = "inner", by = c("Gene1", "Gene2"))
          }
        }
        # Edges at this point may include motif-to-gene links. We need gene-to-gene links.
        if(tf_activity_type %in% c("both", "motif")){
          genegene = DF %>% subset(!(Gene1 %in% names(motif_info$all_motifs)))
          motifgene = DF %>% subset(Gene1 %in% names(motif_info$all_motifs))
          motifgene[["motif"]] = motifgene[["Gene1"]]
          motifgene[["Gene1"]] = NULL
          motifgene %<>% merge(link_genes_to_motifs(motif_info), by = "motif")
          DF = rbind(motifgene[colnames(genegene)], genegene)
        }
        DF$q = DF$knockoff_stat %>% rlookc::knockoffQvals(offset = 0)
        
        cat("\n", as.character(Sys.time()), "\n")
        cat("\nLoading and merging ChIP data...\n")
        mouse_chip = load_chip_data()
        mouse_chip[["is_verified"]] = T
        # This merge makes true edges correct; the rest is NA and fails to distinguish between
        # unknowns and known negatives.
        DF$Gene1 %<>% toupper()
        DF$Gene2 %<>% toupper()
        DF %<>% merge(mouse_chip, all.x = T, all.y = F, by = c("Gene1", "Gene2"))
        # This treats every NA as a known negative -- also wrong but see next steps.
        DF[["is_verified"]][is.na(DF[["is_verified"]])] = FALSE
        # This indicates which hypotheses we could test via ChIP.
        DF[["is_testable"]] = DF[["Gene1"]] %in% mouse_chip[["Gene1"]]
        # This fixes unknown edges (those lacking ChIP data).
        DF[["is_verified"]][ !( DF[["is_testable"]] ) ] = NA
        write.csv(DF, "all_hypotheses.csv")
        
        
        
        # Compute calibration!
        DF = subset(DF, is_testable)
        cat("\n", as.character(Sys.time()), "\n")
        cat("\nChecking calibration \n")
        calibration = list()
        for(fdr in c(1:100)/100){
          calibration[[100*fdr]] =
            DF %>%
            subset(q<fdr) %>%
            dplyr::summarize(
              empirical_fdr = 1-mean(is_verified, na.rm = T),
              num_discoveries = sum(!is.na(is_verified))
            ) %>%
            subset(!is.na(empirical_fdr)) %>%
            dplyr::mutate(nominal_fdr = fdr)
        }
        calibration %<>% data.table::rbindlist()
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
        cat("\n", as.character(Sys.time()), "\n")
        cat("Done.\n")
      })
    })
  })
}




cat("\n", as.character(Sys.time()), "\n")
cat("Starting threads for different experiments. See logs/ for progress.\n")
# Try to re-use any existing results
results = parallel::mclapply(seq(nrow(conditions)), do_one, reuse_results = T, spread_load = T, mc.cores = parallel::detectCores()-1)
print(results)
# Redo failed runs
results2 = parallel::mclapply(which(!sapply(results, is.null)), do_one, mc.cores = parallel::detectCores()-1)
print(results2)
results3 = parallel::mclapply(which(!sapply(results2, is.null)), do_one, mc.cores = parallel::detectCores()-1)
print(results3)
cat("\n", as.character(Sys.time()), "\n")
cat("Done.\n")
