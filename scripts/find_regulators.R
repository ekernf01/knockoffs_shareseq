# Deploy knockoff-based variable selection
# on the TF-target inference problem.
source("../scripts/setup.R")
conditions = read.table(
  header = T,
  text =
    "knockoff_type tf_activity condition_on cell_count_cutoff error_mode seed keratinocyte_only require_motif_support
gaussian rna none 10 none 1 F F
gaussian rna none 100 none 1 F F
gaussian rna none 500 none 1 F F
naive rna none  10 none 1 F F
naive rna none 100 none 1 F F
naive rna none 500 none 1 F F
gaussian rna none 500 downsample 1 F F
gaussian rna none  10 resample 1 F F
gaussian rna none 100 resample 1 F F
gaussian rna none 500 resample 1 F F
gaussian motif none 10 none 1 F F
gaussian motif none 100 none 1 F F
gaussian motif none 500 none 1 F F
")
conditions = rbind(
  conditions, 
  dplyr::mutate(conditions, keratinocyte_only=T)
)


write.csv(conditions, "experiments_to_run.csv")
old_wd = getwd()
condition_idx = 1 # For interactive debugging

do_one = function(condition_idx){
  set.seed(conditions[condition_idx,"seed"])
  Sys.sleep(condition_idx * 120) # This spreads out the peak memory load
  dir.create("logs")
  withr::with_output_sink(file.path("logs", condition_idx), {
    # Set up environment
    attach(conditions[condition_idx,], warn.conflicts = F)
    cat("\nRunning condition:\n")
    write.table(conditions[condition_idx,], quote = F, sep = "\t")
    new_wd = create_experiment_path(conditions, condition_idx)
    # Load data
    normalized_data = set_up_share_skin_pseudobulk(conditions, condition_idx)
    # go to output
    dir.create(new_wd, recursive = T, showWarnings = F)
    withr::with_dir( new_wd, {
      # For this analysis, we want expression centered and scaled
      # For genes with zero variance, esp after downsampling,
      # there's too much brittle software downstream. (Not all of it mine!)
      # So we randomize instead of keeping constants or zeroes.
      cat("\nScaling RNA data...\n")
      safe_scale = function(x, ...){
        if(sd(x)==0){
          return( x + rnorm(length(x)) )
        } else {
          return( (x - mean(x)) / sd(x) )
        }
      }
      normalized_data$mouse_non_tf_expression   %<>% apply(2, safe_scale)
      normalized_data$mouse_tf_expression       %<>% apply(2, safe_scale)
      normalized_data$mouse_tf_expression_noisy %<>% apply(2, safe_scale)
      normalized_data$motif_activity %<>% apply(2, safe_scale)
      # Choose a measure of TF activity
      if(tf_activity == "motif"){
        tf_activity = normalized_data$motif_activity
      } else if(tf_activity=="rna"){
        tf_activity = normalized_data$mouse_tf_expression
      } else{
        stop("tf_activity must be 'rna' or 'motif'")
      }
      # Add simulated measurement error to TF activity
      if(error_mode == "none"){
        tf_activity_before_noise = tf_activity
      } else {
        stopifnot("If tf_activity is 'motif', then error_mode must be 'none'."=tf_activity=="rna")
        tf_activity = normalized_data$mouse_tf_expression_noisy
        tf_activity_before_noise = normalized_data$mouse_tf_expression
      } 
      # Add covariates / surrogate variables
      if(condition_on == "pca"){
        atac_pca = irlba::prcomp_irlba(t(normalized_data$pseudo_bulk_atac), n = 10)
        rna_pca = irlba::prcomp_irlba(normalized_data$pseudo_bulk_rna, n = 10)
        covariates = cbind( tf_activity, atac_pca$x, rna_pca$x)
      } else if(condition_on == "none"){
        covariates = tf_activity
      } else {
        stop("condition_on must be 'pca' or 'none'.")
      }
      # Generate knockoffs
      cat("\nGenerating knockoffs\n")
      if(knockoff_type == "naive"){
        # Assumes variables are in columns! Otherwise you're
        # scrambling all genes within a sample, not all samples within a gene.
        knockoffs = apply(covariates, 2, sample, replace = F)
      } else if (knockoff_type=="gaussian"){
        shrink_by = corpcor::estimate.lambda(covariates)
        Sigma = corpcor::cov.shrink(x = covariates, lambda = shrink_by, lambda.var = 0 )
        Sigma = as(Sigma, "matrix")
        t1 = Sys.time()
        knockoffs = rlookc::createHighDimensionalKnockoffs(
          covariates, lambda = shrink_by,
        )
        t2 = Sys.time()
        t2 - t1
      } else {
        stop("Unknown type of knockoffs requested by upstream code\n.")
      }
      # Strip out extra covariates
      knockoffs = knockoffs[,1:ncol(tf_activity)]
      # Save knockoff realization
      dir.create("output_knockoffs", recursive = T, showWarnings = F)
      saveRDS(knockoffs, "output_knockoffs/knockoffs.Rda")
      saveRDS(tf_activity, "output_knockoffs/original_features.Rda")

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

      # Check calibration with simulated targets
      cat("\nChecking calibration with simulated targets\n")
      dir.create("calibration", showWarnings = F, recursive = T)
      calibration_typical = rlookc::simulateY(
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
      cat("\nChecking calibration with real targets\n")
      cat("Computing knockoff statistics. This may take a while.\n")
      w = lapply(
        seq(ncol(normalized_data$mouse_non_tf_expression)),
        function(i){
          if(i %% 100 == 0){ cat( "\n" ); cat(i); cat(" ")}
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
      DF_tf_target = list()
      for(k in seq(ncol(normalized_data$mouse_non_tf_expression))){
        DF_tf_target[[k]] = data.frame(
          Gene1 = colnames(normalized_data$mouse_tf_expression),
          Gene2 = colnames(normalized_data$mouse_non_tf_expression)[ k],
          knockoff_stat = w[[k]]
        ) %>%
          subset(abs(knockoff_stat) > 0)
      }
      cat("Assembling knockoff stats and cleaning up.\n")
      rm("w"); gc()
      DF_tf_target = data.table::rbindlist(DF_tf_target)
      # Save results
      write.csv(DF_tf_target, "output_knockoffs/knockoff_stats.csv", append = T, col.names = F)
      rm("DF_tf_target"); gc()

      # We will check FDR for a subset of TF-target hypotheses supported by a motif in a region whose chromatin accessibility
      # correlates with the target gene expression.
      gene_coords = biomaRt::getBM(attributes = c("mgi_symbol", "chromosome_name", "start_position", "end_position"),
                          filters = "mgi_symbol",
                          values = normalized_data$gene_metadata$Gene1,
                          mart = biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl"))
      gene_coords = GRanges(
        seqnames = paste0("chr", gene_coords$chromosome_name),
        ranges = IRanges(gene_coords$start_position, gene_coords$end_position),
        genome = "mm10"
      )
      skin_atac_peaks = GRanges(
        seqnames = skin_atac_peaks$V1,
        ranges = IRanges(skin_atac_peaks$V2, skin_atac_peaks$V3),
        genome = "mm10"
      )
      overlaps = GenomicRanges::findOverlaps(
        gene_coords,
        skin_atac_peaks,
        maxgap = 1e5
      )
      overlaps %<>% as.data.frame() %>% set_colnames(c("gene", "enhancer"))
      overlaps[["distance"]] = GenomicRanges::distance(
        gene_coords[overlaps[["gene"]]],
        skin_atac_peaks[overlaps[["enhancer"]]]
      )
      overlaps$correlation = NA
      for(i in seq_along(overlaps[[1]])){
        if((i%%10000)==0){cat(".")}
        peak_idx = overlaps[i, "enhancer"]
        gene_idx = overlaps[i, "gene"]
        overlaps[i, "correlation"] = cor(
            normalized_data$pseudo_bulk_atac[peak_idx,],
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
      motif_enhancer_links = summary(motif_info$motif_ix@assays@data@listData$motifMatches) %>% set_colnames(c("enhancer", "motif_index", "is_connected"))
      gene_motif_links = motif_info$all_motifs %>%
        sapply(TFBSTools::name) %>%
        data.frame(genes=.) %>%
        (dplyr::add_rownames) %>%
        dplyr::rename(motif_index = rowname) %>%
        dplyr::mutate(genes = gsub("\\(var\\..*\\)", "", genes)) %>%
        tidyr::separate_rows("genes", sep = "::") %>%
        dplyr::mutate(Gene1 = stringr::str_to_sentence(genes))
      motif_gene_links = merge(motif_enhancer_links, enhancer_gene_links, by = "enhancer")
      gene_gene_links_from_motif_analysis = merge(gene_motif_links, motif_gene_links, by = "motif_index")
      gene_gene_links_from_motif_analysis = gene_gene_links[c("Gene1", "Gene2")]

      # Load chip-atlas target gene lists
      cat("\nLoading ChIP data...\n")
      withr::with_dir(
        DATALAKE,
        {
          average_signals = function(DF){
            new_df = DF[1:2]
            new_df$average_signal = DF[-(1:2)] %>% rowMeans()
            return(new_df)
          }
          chip_files = list.files("chip-atlas/filtered_by_celltype/skin/hg19", full = T)[1:3]
          mouse_chip =
            lapply(chip_files, read.csv, sep = " ", header = T) %>%
            lapply(average_signals) %>%
            lapply(subset, average_signal>mean(average_signal)) %>% # filter for strong signals
            data.table::rbindlist() %>%
            dplyr::mutate(is_verified = T)
        }
      )

      # Load data and look only at the motif-supported hypotheses
      DF = read.csv("output_knockoffs/knockoff_stats.csv")
      # This merge makes true edges correct; the rest is NA and fails to distinguish between
      # unknowns and known negatives.
      DF %<>% merge(mouse_chip, all.x = T, all.y = F)
      # This treats every NA as a known negative -- also wrong but see next steps.
      DF[["is_verified"]][is.na(DF[["is_verified"]])] = FALSE
      # This indicates which hypotheses we could test via ChIP.
      DF[["is_testable"]] = DF[["Gene1"]] %in% mouse_chip[["Gene1"]]
      # This fixes unknown edges (those lacking ChIP data).
      DF[["is_verified"]][ !( DF[["is_testable"]] ) ] = NA

      # Check FDR on 1) testable hypotheses, with 2) motif support
      DF = subset(DF, is_testable)
      if( require_motif_support ){
        DF = merge(DF, gene_gene_links_from_motif_analysis, type = "inner")
      }
      DF$q = DF$knockoff_stat %>% rlookc::knockoffQvals(offset = 0)

      # Compute calibration!
      calibration = list()
      for(fdr in c(1:100)/100){
        calibration[[100*fdr]] =
          DF %>%
          subset(q<fdr) %>%
          dplyr::group_by(Gene1) %>%
          dplyr::summarize(
            empirical_fdr = 1-mean(is_verified, na.rm = T),
            num_discoveries = sum(!is.na(is_verified))
          ) %>%
          subset(!is.na(empirical_fdr)) %>%
          dplyr::mutate(nominal_fdr = fdr)
      }
      calibration %<>% data.table::rbindlist()
      # Add simple binomial standard errors
      calibration %<>% dplyr::mutate( moe_95 = 1.96 * sqrt( empirical_fdr * ( 1 - empirical_fdr ) / num_discoveries ) )
      write.csv(calibration, "calibration/chip_calibration.csv")

      # Plot and save calibration
      ggplot(calibration) +
        geom_point(aes(x = nominal_fdr, y = empirical_fdr, colour = Gene1, shape = Gene1)) +
        geom_errorbar(aes(x = nominal_fdr, colour = Gene1,
                          ymin = pmax(0, empirical_fdr - moe_95),
                          ymax = pmin(1, empirical_fdr + moe_95)
        )) +
        ggtitle("Calibration vs ChIP-Atlas") +
        scale_y_continuous(limits = 0:1)
      ggsave("calibration/chip_calibration.pdf", width = 4, height = 4)
      cat("Done.\n")
    })
  })
}
cat("Starting threads for different experiments. See logs/ for progress.\n")
results = parallel::mclapply(seq(nrow(conditions)), do_one, mc.cores = parallel::detectCores()-1)
saveRDS(results, "logs/mclapply_output.Rdata")
conditions %>%
  cbind(successfully_ran = sapply(results, is.null)) %>%
  write.table(sep = "\t", quote = F)
# Redo failed runs
parallel::mclapply(which(!sapply(results, is.null)), do_one, mc.cores = parallel::detectCores()-1)
cat("Done.\n")
