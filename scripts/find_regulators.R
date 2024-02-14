# Deploy knockoff-based variable selection
# on the TF-target inference problem.
source("../scripts/setup.R")
conditions = read.table(
  header = T,
  text =
    "knockoff_type tf_activity_type condition_on cell_count_cutoff error_mode seed celltype require_motif_support only_motif_support include_decoys
gaussian   rna none  10       none 1 skin F  F  F
gaussian   rna none 100       none 1 skin F  F  F
gaussian   rna none 500       none 1 skin F  F  F
permuted   rna none  10       none 1 skin F  F  F
permuted   rna none 100       none 1 skin F  F  F
permuted   rna none 500       none 1 skin F  F  F
gaussian   rna none  10   resample 1 skin F  F  F
gaussian   rna none 100   resample 1 skin F  F  F
gaussian   rna none 500   resample 1 skin F  F  F
gaussian motif none  10       none 1 skin F  F  F
gaussian motif none 100       none 1 skin F  F  F
gaussian motif none 500       none 1 skin F  F  F
gaussian  both none  10       none 1 skin F  F  F
gaussian  both none 100       none 1 skin F  F  F
gaussian  both none 500       none 1 skin F  F  F
gaussian  both  pca  10       none 1 skin F  F  F
gaussian  both  pca 100       none 1 skin F  F  F
gaussian  both  pca 500       none 1 skin F  F  F
gaussian  both  pca  10       none 1 skin  T F  F
gaussian  both  pca 100       none 1 skin  T F  F
gaussian  both  pca 500       none 1 skin  T F  F
gaussian  both  pca  10       none 1 skin F   T F
gaussian  both  pca 100       none 1 skin F   T F
gaussian  both  pca 500       none 1 skin F   T F
")
conditions = Reduce(
  rbind, 
  list(
    conditions, 
    dplyr::mutate(conditions, celltype="keratinocyte"),
    dplyr::mutate(subset(conditions, cell_count_cutoff < 500), celltype="pbmc"),
    dplyr::mutate(subset(conditions, cell_count_cutoff < 500), celltype="tcell"),
    dplyr::mutate(subset(conditions, cell_count_cutoff < 500), celltype="pbmc_subset")
  )
)
rownames(conditions) = NULL

conditions %<>% rbind(
  read.table(
  header = T,
  text =
    "knockoff_type tf_activity_type condition_on cell_count_cutoff error_mode seed celltype require_motif_support only_motif_support include_decoys
gaussian   rna none  50       none 1 tcell F  F  F
gaussian   rna none 200       none 1 tcell F  F  F
permuted   rna none  50       none 1 tcell F  F  F
permuted   rna none 200       none 1 tcell F  F  F
gaussian   rna none  50   resample 1 tcell F  F  F
gaussian   rna none 200   resample 1 tcell F  F  F
gaussian motif none  50       none 1 tcell F  F  F
gaussian motif none 200       none 1 tcell F  F  F
gaussian  both none  50       none 1 tcell F  F  F
gaussian  both none 200       none 1 tcell F  F  F
"), .
)
write.csv(conditions, "experiments_to_run.csv")
old_wd = getwd()
condition_idx = 1 # For interactive debugging

#' Run a single experiment described by a row of the 'conditions' dataframe defined above. Its columns are documented here as if they were parameters.
#'
#' @param knockoff_type 'permuted' or 'gaussian'. 
#' @param tf_activity 'rna' or 'motif' or 'both'. How to extract TF activities. 
#' @param condition_on 'none' or 'pca.' if 'pca', then include the top 10 PC's from both ATAC and RNA data.
#' @param cell_count_cutoff If a cluster has fewer cells than this, it gets excluded from our analysis.
#' @param error_mode If 'none', take no special action. If 'resample', replace mRNA count Xij with a Poisson(Xij). 
#' @param seed Random seed in case you want to replicate experiments. 
#' @param celltype skin: all the mouse skin SHARE data. keratinocyte: a subset of the mouse skin SHARE data. pbmc: the human data. tcell: a subset of the pbmc data.
#' @param require_motif_support If T, retain an edge only when there is a motif for the TF in a nearby open chromatin region.
#' @param only_motif_support If T, ignore everything else and just get the FDR you'd see if you used motif-matching as the sole inference method.
#' 
do_one = function(condition_idx, reuse_results = F){
  set.seed(conditions[condition_idx,"seed"])
  dir.create("logs", showWarnings = F)
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
        if(!reuse_results){
          normalized_data$tf_expression %>% dim %>% write.csv("tf_expression_dimensions.csv")
          
          # For this analysis, we want expression centered and scaled
          # For genes with zero variance, esp after downsampling,
          # there's too much brittle software downstream. (Not all of it mine!)
          # So we randomize instead of keeping constants or zeroes.
          cat("\n", as.character(Sys.time()), "\n")
          cat("\nScaling RNA data...\n")
          
          normalized_data$non_tf_expression   %<>% apply(2, safe_scale)
          normalized_data$tf_expression       %<>% apply(2, safe_scale)
          normalized_data$tf_expression_noisy %<>% apply(2, safe_scale)
          normalized_data$motif_activity %<>% apply(2, safe_scale)
          # Add decoy TF's: correlated at r^2 = 0.9 with 100 randomly chosen TF's, but conditionally independent from targets. 
          # R^2=0.9 is determined by the factor of 1/3 in the code.
          if(include_decoys){
            decoy_parents = sample(colnames(normalized_data$tf_expression), 100, replace = F)
            decoys = normalized_data$tf_expression[, decoy_parents] 
            decoys = decoys + (1/3)*matrix( rnorm(prod(dim(decoys))), nrow = nrow(decoys), ncol = ncol(decoys) )
            colnames(decoys) = paste0("DECOY_", decoy_parents)
            normalized_data$tf_expression %<>% cbind(decoys)
            normalized_data$tf_expression_noisy %<>% cbind(decoys)
          }
          
          # Choose a measure of TF activity
          if(tf_activity_type == "motif"){
            tf_activity = normalized_data$motif_activity
            colnames(tf_activity) = colnames(normalized_data$motif_activity)
            differential_motifs = 
              normalized_data$motif_activity %>% 
              scale %>% 
              data.frame %>%
              tibble::rownames_to_column("cluster") %>%
              tidyr::pivot_longer(matches("MA"), names_to = "motif", values_to = "Zscore") %>%
              dplyr::group_by(cluster) %>%
              dplyr::top_n(Zscore, n=5) %>%
              dplyr::ungroup() %>%
              dplyr::mutate(cluster = gsub("X", "", cluster)) %>%
              merge(link_genes_to_motifs(motif_info)) %>% 
              merge(normalized_data$pseudo_bulk_metadata) %>%
              dplyr::arrange(largest_subpopulation, -Zscore) 
            differential_motifs %>%
              write.csv("motif_score_enrichment.csv")
            
          } else if(tf_activity_type=="both"){
            tf_activity =       cbind(         normalized_data$motif_activity,           normalized_data$tf_expression)
            colnames(tf_activity) = c(colnames(normalized_data$motif_activity), colnames(normalized_data$tf_expression))
          } else if(tf_activity_type=="rna"){
            tf_activity =                    normalized_data$tf_expression
            colnames(tf_activity) = colnames(normalized_data$tf_expression)
          } else{
            stop("tf_activity_type must be 'rna' or 'motif' or 'both'")
          }
          # Add simulated measurement error to TF activity
          if(error_mode == "none"){
            tf_activity_before_noise = tf_activity
          } else {
            stopifnot("Unless tf_activity_type is 'rna', error_mode must be 'none'."=tf_activity_type=="rna")
            tf_activity              = normalized_data$tf_expression_noisy
            tf_activity_before_noise = normalized_data$tf_expression
          } 
          # Add covariates / surrogate variables
          if(condition_on == "pca"){
            atac_pca = irlba::prcomp_irlba(normalized_data$pseudo_bulk_atac, n = 10)
            rna_pca  = irlba::prcomp_irlba(normalized_data$pseudo_bulk_rna, n = 10)
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
          w = parallel::mclapply(
            seq(ncol(normalized_data$non_tf_expression)),
            function(i){
              if(i %% 100 == 0){ cat( "\n", i, "/", ncol(normalized_data$non_tf_expression) )}
              fast_lasso_penalty(
                y   = normalized_data$non_tf_expression[,i],
                X   = tf_activity,
                X_k = knockoffs
              )
            }, 
            mc.cores = parallel::detectCores()-1, 
            mc.preschedule = F
          )
          length(w[[1]])
          length(w)
          cat("Saving knockoff statistics.\n")
          saveRDS(w, "output_knockoffs/w.Rds")
          
          # Assemble results nicely
          cat("\n", as.character(Sys.time()), "\n")
          cat("Assembling knockoff stats and cleaning up.\n")
          DF_tf_target = list()
          for(k in seq(ncol(normalized_data$non_tf_expression))){
            if(k %% 100 == 0){ cat( "\n", k, "/", ncol(normalized_data$non_tf_expression) )}
            DF_tf_target[[k]] = data.frame(
              Gene1 = colnames(tf_activity),
              Gene2 = colnames(normalized_data$non_tf_expression)[ k],
              knockoff_stat = w[[k]]
            ) %>%
              subset(abs(knockoff_stat) > 0)
          }
          rm("w"); gc()
          DF_tf_target = data.table::rbindlist(DF_tf_target)
          
          # Save results
          write.csv(DF_tf_target, "output_knockoffs/knockoff_stats.csv")
          rm("DF_tf_target"); gc()
          
          
          # We will try using 1) all testable hypotheses, with 2) testable hypotheses motif support, or 3) as a control, only motif-based hypotheses
          cat("\n", as.character(Sys.time()), "\n")
          cat("\nIntegrating knockoff-based and motif-based hypotheses... \n")
          if( only_motif_support ){
            DF = get_motif_supported_hypotheses(normalized_data, celltype = celltype)
            DF$Gene1 %<>% toupper()
            DF$Gene2 %<>% toupper()
            DF$knockoff_stat = 100 + rnorm(nrow(DF)) #these fake knockoff statistics will yield q=0 for all hypotheses
          } else {
            DF = read.csv("output_knockoffs/knockoff_stats.csv")
            DF$Gene1 %<>% toupper()
            DF$Gene2 %<>% toupper()
            if( require_motif_support ){
              mosh = get_motif_supported_hypotheses(normalized_data, celltype = celltype)
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
          all_chip = load_chip_data(celltype)
          all_chip[["is_verified"]] = T
          # This merge makes true edges correct; the rest is NA and fails to distinguish between
          # unknowns and known negatives.
          DF$Gene1 %<>% toupper()
          DF$Gene2 %<>% toupper()
          DF %<>% merge(all_chip, all.x = T, all.y = F, by = c("Gene1", "Gene2"))
          # This treats every NA as a known negative -- also wrong but see next steps.
          DF[["is_verified"]][is.na(DF[["is_verified"]])] = FALSE
          # This indicates which hypotheses we could test via ChIP.
          DF[["is_testable"]] = DF[["Gene1"]] %in% all_chip[["Gene1"]]
          # This fixes unknown edges (those lacking ChIP data).
          DF[["is_verified"]][ !( DF[["is_testable"]] ) ] = NA
          write.csv(DF, gzfile("all_hypotheses.csv.gz"))
        }
        
        # Compute calibration
        DF = read.csv("all_hypotheses.csv.gz")
        DF = subset(DF, is_testable)
        cat("\n", as.character(Sys.time()), "\n")
        cat("\nChecking calibration vs chip \n")
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
        dir.create("calibration", showWarnings = F)
        write.csv(calibration, "calibration/chip_calibration.csv")
        
        DF = read.csv("all_hypotheses.csv.gz")
        cat("\n", as.character(Sys.time()), "\n")
        cat("\nChecking calibration via decoys \n")
        calibration = list()
        for(fdr in c(1:100)/100){
          calibration[[100*fdr]] =
            DF %>%
            subset( q < fdr ) %>%
            dplyr::mutate( is_decoy = grepl( "^DECOY", Gene1 ) ) %>%
            dplyr::summarize(
              empirical_fdr = mean( is_decoy, na.rm = T ),
              num_discoveries = sum( !is.na( is_decoy ) )
            ) %>%
            subset(!is.na(empirical_fdr)) %>%
            dplyr::mutate(nominal_fdr = fdr)
        }
        calibration %<>% data.table::rbindlist()
        calibration %<>% dplyr::mutate( moe_95 = 1.96 * sqrt( empirical_fdr * ( 1 - empirical_fdr ) / num_discoveries ) )
        dir.create("calibration", showWarnings = F)
        write.csv(calibration, "calibration/decoy_calibration.csv")
        
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
  return()
}

cat("\n", as.character(Sys.time()), "\n")
# Try to re-use any existing results (this is useful for interacting with the code)
results = lapply(
  seq(nrow(conditions)),
  do_one, 
  reuse_results = F)
print(results)
# Redo failed runs
todo = which(!sapply(results, is.null))
results = lapply(todo, do_one)
print(results)
cat("\n", as.character(Sys.time()), "\n")
cat("Done.\n")
