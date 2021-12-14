# Deploy certain variable selection techniques with unique advantages
# on the TF-target inference problem.
source("../scripts/setup.R")
conditions = read.table(
  header = T,
  text =
    "knockoff_type condition_on cell_count_cutoff error_mode seed keratinocyte_only
naive none  10 none 1 F
naive none 100 none 1 F
naive none 500 none 1 F
gaussian none  10 none 1 F
gaussian none 100 none 1 F
gaussian none 500 none 1 F
gaussian none  10 resample 1 F
gaussian none 100 resample 1 F
gaussian none 500 resample 1 F
gaussian none 500 downsample 1 F
gaussian atac  10 none 1 F
gaussian atac 100 none 1 F
gaussian atac 500 none 1 F
gaussian rna5  10 none 1 F
gaussian rna5 100 none 1 F
gaussian rna5 500 none 1 F
gaussian rna10  10 none 1 F
gaussian rna10 100 none 1 F
gaussian rna10 500 none 1 F
gaussian rna15  10 none 1 F
gaussian rna15 100 none 1 F
gaussian rna15 500 none 1 F
gaussian atac  10 none 1 T
gaussian atac 100 none 1 T
gaussian atac 500 none 1 T
gaussian rna5  10 none 1 T
gaussian rna5 100 none 1 T
gaussian rna5 500 none 1 T
gaussian rna10  10 none 1 T
gaussian rna10 100 none 1 T
gaussian rna10 500 none 1 T
gaussian rna15  10 none 1 T
gaussian rna15 100 none 1 T
gaussian rna15 500 none 1 T"
)
write.csv(conditions, "experiments_to_run.csv")
old_wd = getwd()

do_one = function(condition_idx){
  set.seed(conditions[condition_idx,"seed"])
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
      # there's too much brittle software downstream. (Not all mine! Heehee!)
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

      # Get ATAC PC's
      if(condition_on == "atac"){
        atac_pca = irlba::prcomp_irlba(t(normalized_data$pseudo_bulk_atac), n = 10)
        covariates = cbind( normalized_data$mouse_tf_expression_noisy, atac_pca$x)
      } else if(condition_on == "rna5"){
        rna_pca = irlba::prcomp_irlba(normalized_data$pseudo_bulk_rna, n = 5)
        covariates = cbind( normalized_data$mouse_tf_expression_noisy, rna_pca$x)
      } else if(condition_on == "rna10"){
        rna_pca = irlba::prcomp_irlba(normalized_data$pseudo_bulk_rna, n = 10)
        covariates = cbind( normalized_data$mouse_tf_expression_noisy, rna_pca$x)
      } else if(condition_on == "rna15"){
        rna_pca = irlba::prcomp_irlba(normalized_data$pseudo_bulk_rna, n = 15)
        covariates = cbind( normalized_data$mouse_tf_expression_noisy, rna_pca$x)
      } else {
        covariates = normalized_data$mouse_tf_expression_noisy
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
      knockoffs = knockoffs[,1:ncol(normalized_data$mouse_tf_expression)]
      # Save knockoff realization
      dir.create("output_knockoffs", recursive = T, showWarnings = F)
      saveRDS(knockoffs, "output_knockoffs/knockoffs.Rda")
      saveRDS(normalized_data$mouse_tf_expression_noisy, "output_knockoffs/original_features.Rda")

      # This is a fast way to compute variable importance statstics.
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

      # Check calibration
      cat("\nChecking calibration with simulated targets\n")
      dir.create("calibration", showWarnings = F, recursive = T)
      calibration_typical = rlookc::simulateY(
        n_sim = 100,
        X = normalized_data$mouse_tf_expression,
        knockoffs = knockoffs,
        statistic = fast_lasso_penalty,
        plot_savepath = "calibration/average_case_calibration_X_exact.pdf"
      )
      write.csv(calibration_typical$calibration$fdr, "calibration/average_case_calibration_X_exact.csv")

      if(error_mode != "none"){
        cat("\nChecking calibration with simulated targets and specified errors in X \n")
        dir.create("calibration", showWarnings = F, recursive = T)
        calibration_typical = rlookc::simulateY(
          n_sim = 100,
          X = normalized_data$mouse_tf_expression,
          X_observed = normalized_data$mouse_tf_expression_noisy,
          knockoffs = knockoffs,
          statistic = fast_lasso_penalty,
          plot_savepath = "calibration/average_case_calibration_X_errors.pdf"
        )
        write.csv(calibration_typical$calibration$fdr, "calibration/average_case_calibration_X_errors.csv")

      }
      # Find regulators of all genes.
      # Deploy on non-TF targets
      cat("\nChecking calibration with real targets\n")
      cat(
        "Length check:",
        length(
          fast_lasso_penalty(
            y   = normalized_data$mouse_non_tf_expression[,1],
            X   = normalized_data$mouse_tf_expression,
            X_k = knockoffs)
        )
      )
      cat("Computing knockoff statistics. This may take a while.\n")
      w = lapply(
        seq(ncol(normalized_data$mouse_non_tf_expression)),
        function(i){
          if(i %% 100 == 0){ cat( "\n" ); cat(i); cat(" ")}
          fast_lasso_penalty(
            y   = normalized_data$mouse_non_tf_expression[,i],
            X   = normalized_data$mouse_tf_expression,
            X_k = knockoffs
          )
        }
      )
      length(w[[1]])
      length(w)
      # Deploy on TF targets
      w_tf = list()
      for(i in seq(ncol(normalized_data$mouse_tf_expression))){
        if(i %% 100 == 0){ cat( "\n" ); cat(i); cat(" ")}
        # knockoffs_minus_i = rlookc::formOneLook(
        #   knockoffs = looks_compact$knockoffs,
        #   updates = looks_compact$updates,
        #   vars_to_omit = looks_compact$vars_to_omit,
        #   k = i
        # )
        knockoffs_minus_i = knockoffs[,-i]
        w_tf[[i]] = fast_lasso_penalty(
          y   = normalized_data$mouse_tf_expression[,i],
          X   = normalized_data$mouse_tf_expression[,-i],
          X_k = knockoffs_minus_i
        )
      }
      cat("Saving knockoff statistics.\n")
      saveRDS(w, "output_knockoffs/w.Rds")
      saveRDS(w_tf, "output_knockoffs/wtf.Rds")

      # Assemble results nicely
      DF_tf_target = DF_tf_tf = list()
      for(k in seq(ncol(normalized_data$mouse_tf_expression))){
        DF_tf_tf[[k]] = data.frame(
          Gene1 = colnames(normalized_data$mouse_tf_expression)[-k],
          Gene2 = colnames(normalized_data$mouse_tf_expression)[ k],
          knockoff_stat = w_tf[[k]]
        ) %>%
          subset(abs(knockoff_stat) > 0)
      }
      for(k in seq(ncol(normalized_data$mouse_non_tf_expression))){
        DF_tf_target[[k]] = data.frame(
          Gene1 = colnames(normalized_data$mouse_tf_expression),
          Gene2 = colnames(normalized_data$mouse_non_tf_expression)[ k],
          knockoff_stat = w[[k]]
        ) %>%
          subset(abs(knockoff_stat) > 0)
      }
      cat("Assembling knockoff stats and cleaning up.\n")
      rm("w", "w_tf"); gc()
      DF_tf_tf = data.table::rbindlist(DF_tf_tf)
      DF_tf_target = data.table::rbindlist(DF_tf_target)
      write.csv(DF_tf_tf, "output_knockoffs/knockoff_stats.csv")
      rm("DF_tf_tf"); gc()
      write.csv(DF_tf_target, "output_knockoffs/knockoff_stats.csv", append = T, col.names = F)
      rm("DF_tf_target"); gc()
      DF = read.csv("output_knockoffs/knockoff_stats.csv")
      DF$q = DF$knockoff_stat %>% rlookc::knockoffQvals(offset = 0)
      # How much signal is there?
      {
        pdf("q_value_distribution.pdf")
        plot(ecdf(DF$q), xlim = 0:1)
        hist(DF$q, xlab = "Knockoff q-values", ylab = "How many hits", main = "\"Power\"")
        dev.off()
      }

      write.csv(
        file = "q_value_distribution.csv",
        x = data.frame(
          fdr_cutoff = (1:5)/10,
          n_discoveries = ((1:5)/10) %>% sapply(function(fdr) sum(DF$q<fdr)
          )
        )
      )

      # Load chip-atlas target gene lists
      cat("\nLoading ChIP data...\n")
      withr::with_dir(
        DATALAKE,
        {
          chip_meta = read.table("chip-atlas/mouse/experimentList.tab.standard.fields.only.tab",
                                 sep = "\t", header = F)
          colnames(chip_meta) =  c("Experimental ID",
                                   "Genome assembly"	,
                                   "Antigen class",
                                   "Antigen",
                                   "Cell type class",
                                   "Cell type",
                                   "Cell type description",
                                   "Processing logs")
          chip_meta = subset(chip_meta, `Genome assembly`=="mm10")
          chip_meta %<>% subset(`Antigen class`=="TFs and others")
          chip_meta %<>% subset(`Cell type class`=="Epidermis")
          chip_files =
            list.files("chip-atlas/mouse/targets", full = T)
          names(chip_files) = gsub(".10.tsv$", "", basename(chip_files))
          chip_files = chip_files[unique(chip_meta[["Antigen"]])]
          mouse_chip =
            lapply(chip_files, read.csv, sep = "\t", header = T) %>%
            lapply(extract2, 1) %>%
            mapply(
              function(X, tf) {
                data.frame("Gene1" = gsub(".10.tsv$", "", basename(tf)), "Gene2" = X)
              },
              X = .,
              tf = chip_files,
              SIMPLIFY = F) %>%
            data.table::rbindlist() %>%
            dplyr::mutate(is_verified = T)
        }
      )

      # This makes true edges correct; the rest is NA and fails to distinguish between
      # unknowns and known negatives.
      DF %<>% merge(mouse_chip, all.x = T, all.y = F)
      # This treats every NA as a known negative.
      DF[["is_verified"]][is.na(DF[["is_verified"]])] = F
      # This fixes unknown edges (those lacking ChIP data).
      DF[["is_verified"]][ !( DF[["Gene1"]] %in% mouse_chip[["Gene1"]] ) ] = NA


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
results = parallel::mclapply(seq(nrow(conditions)), do_one, mc.cores = parallel::detectCores())
saveRDS(results, "logs/mclapply_output.Rdata")
conditions %>%
  cbind(successfully_ran = sapply(results, is.null)) %>%
  write.table(sep = "\t", quote = F)
cat("Done.\n")
