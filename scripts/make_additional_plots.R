# setwd("~/Desktop/jhu/research/projects/knockoffs/applications/share-seq/v12")
source("../scripts/setup.R")
sanitize_names = function(df){
  colnames(df) = colnames(df) %>% 
    gsub("empirical_fdr", "observed_fdr", .)  %>% 
    gsub("nominal_fdr", "expected_fdr", .) 
  return(df)
}
ggplot2::theme_update(text = element_text(family = "ArialMT"))
conditions = read.csv("experiments_to_run.csv")
# Plot big calibration complicated many facets such oof wow
all_calibration = knn_test = list()
for(condition_idx in seq(nrow(conditions))){
  # Set up condition-specific environment
  attach(conditions[condition_idx,], warn.conflicts = F)
  cat("\nRetrieving condition:\n")
  write.table(conditions[condition_idx,], quote = F, sep = "\t")
  new_wd = create_experiment_path(conditions, condition_idx)
  normalized_data = set_up_share_skin_pseudobulk(conditions, condition_idx)
  # These help condense results from rlookc::simulateY()
  summarize_simulation = function(X) {
    X %>%
      colMeans %>%
      (function(x) data.frame(expected_fdr = gsub("^X", "", names(x)) %>% as.numeric, observed_fdr = x)) %>%
      dplyr::mutate( num_discoveries  = NA,  moe_95 = NA)
  }
  read_simulation_if_exists = function(){
    if(file.exists("calibration/average_case_calibration_X_errors.csv")){
      read.csv("calibration/average_case_calibration_X_errors.csv", row.names = 1) %>%
        summarize_simulation %>%
        cbind(evaluation="simulation_X_has_errors") %>%
        cbind(Gene1="") %>%
        return
    } else {
      return(NULL)
    }
  }
  # Run KNN exchangeability diagnostic
  # Set n_neighbors low because there's not many observations
  withr::with_dir(new_wd, {
    knn_test[[condition_idx]] = rlookc::KNNTest(X = readRDS("output_knockoffs/original_features.Rda"),
                                                X_k = readRDS("output_knockoffs/knockoffs.Rda"),
                                                n_neighbors = 5)
  })
  # Grab outputs and make them uniform
  withr::with_dir(new_wd, {
    all_calibration[[condition_idx]] = list(
      read.csv("calibration/chip_calibration.csv", row.names = 1) %>%
        cbind(evaluation="chip") %>% 
        sanitize_names,
      read.csv("calibration/average_case_calibration_X_exact.csv", row.names = 1) %>%
        summarize_simulation %>%
        sanitize_names %>%
        cbind(evaluation="simulation_X_is_exact") %>%
        cbind(Gene1=""),
      read_simulation_if_exists()
    ) %>%
      (function(X) (X[!sapply(X, is.null)])) %>%
      Reduce(f=rbind) %>%
      merge(conditions[condition_idx,])
  })
}
all_calibration %<>% Reduce(f = rbind)
all_calibration %<>% dplyr::mutate( evaluation = paste(Gene1, evaluation) )

# Plot naive versus "shrinkage" (corpcor gaussian) knockoffs
all_calibration %>%
  subset(error_mode == "none"  & !keratinocyte_only & condition_on =="none" & evaluation == " simulation_X_is_exact") %>%
  ggplot() +
  geom_abline(aes(slope = 1, intercept = 0)) +
  geom_point(aes(x = expected_fdr, y = observed_fdr, color = knockoff_type, shape = knockoff_type)) +
  facet_grid( ~ cell_count_cutoff) +
  ggtitle("Conditional independence testing on SHARE-seq data", "Real TF expression and simulated targets") + 
  theme(text = element_text(family = "ArialMT")) +  
  scale_x_continuous(breaks = (0:2)/2, limits = 0:1) +  
  scale_y_continuous(breaks = (0:2)/2, limits = 0:1)
ggsave("shareseq_naive_vs_gaussian.pdf", width = 5, height = 2)
ggsave("shareseq_naive_vs_gaussian.svg", width = 5, height = 2)

# Plot initial ChIP results
all_calibration %>%
  subset(
    knockoff_type == "gaussian"  &
      !keratinocyte_only &
      condition_on =="none" &
      error_mode == "none" &
      grepl("chip", evaluation) &
      cell_count_cutoff == 10
  ) %>%
  ggplot() +
  geom_abline(aes(slope = 1, intercept = 0)) +
  geom_errorbar(aes(x = expected_fdr, y = observed_fdr, ymin = observed_fdr - moe_95, ymax = observed_fdr + moe_95 ), color = "grey") +
  geom_point(aes(x = expected_fdr, y = observed_fdr, ymin = observed_fdr - moe_95, ymax = observed_fdr + moe_95 )) +
  facet_grid(~evaluation) +
  ggtitle("TRN inference on SHARE-seq data", "Real data") + 
  theme(text = element_text(family = "ArialMT")) +  
  scale_x_continuous(breaks = (0:2)/2, limits = 0:1) +  
  scale_y_continuous(breaks = (0:2)/2, limits = 0:1)
ggsave("shareseq_chip_initial.pdf", width = 5, height = 2)
ggsave("shareseq_chip_initial.svg", width = 5, height = 2)

# Plot results with varied cell count cutoff
all_calibration %>%
  dplyr::mutate(cell_count_cutoff = as.character(cell_count_cutoff)) %>%
  subset(knockoff_type == "gaussian"  & !keratinocyte_only & condition_on =="none" & grepl("chip", evaluation) & error_mode == "none") %>%
  ggplot() +
  geom_abline(aes(slope = 1, intercept = 0)) +
  geom_point(aes(x = expected_fdr, y = observed_fdr, color = cell_count_cutoff, shape = cell_count_cutoff)) +
  facet_grid(~evaluation) +
  ggtitle("Calibration versus minimum required cell count", "Real data") + 
  scale_color_discrete(name = "Cell\ncount\ncutoff") +
  scale_shape_discrete(name = "Cell\ncount\ncutoff") +
  theme(text = element_text(family = "ArialMT")) +  
  scale_x_continuous(breaks = (0:2)/2, limits = 0:1) +  
  scale_y_continuous(breaks = (0:2)/2, limits = 0:1)
ggsave("shareseq_cellcount_cutoff.pdf", width = 5, height = 2)
ggsave("shareseq_cellcount_cutoff.svg", width = 5, height = 2)

# Plot results with/without error in X
all_calibration %>%
  subset(knockoff_type == "gaussian"  &
           !keratinocyte_only &
           condition_on =="none" &
           evaluation != " simulation_X_has_errors" &
           error_mode != "downsample") %>%
  dplyr::mutate(evaluation = gsub(" simulation_X_has_errors", "simulation", evaluation)) %>%
  ggplot() +
  geom_abline(aes(slope = 1, intercept = 0)) +
  geom_point(aes(x = expected_fdr, y = observed_fdr, color = error_mode, shape = error_mode)) +
  facet_grid(cell_count_cutoff ~ evaluation) +
  ggtitle("Calibration with additional Poisson error", "Real and resampled regulators; Real and simulated targets") +
  theme(text = element_text(family = "ArialMT")) +  
  scale_x_continuous(breaks = ((0:2)/2) %>% setNames(c("0", "0.5", "1")), limits = 0:1) +  
  scale_y_continuous(breaks = (0:2)/2, limits = 0:1)
ggsave("shareseq_measurement_error.pdf", width = 6, height = 3.5)
ggsave("shareseq_measurement_error.svg", width = 6, height = 3.5)

# Show results with/without conditioning on confounders
all_calibration %>%
  subset(knockoff_type == "gaussian" & error_mode =="none" & grepl("chip", evaluation) & cell_count_cutoff==10) %>%
  ggplot() +
  geom_abline(aes(slope = 1, intercept = 0)) +
  geom_point(aes(x = expected_fdr, y = observed_fdr, color = condition_on, shape = condition_on)) +
  facet_grid( ifelse(keratinocyte_only, "keratinocyte", "all") ~ evaluation ) +
  ggtitle("Calibration by cell type", "Real data") + 
  theme(text = element_text(family = "ArialMT")) +  
  scale_x_continuous(breaks = ((0:2)/2) %>% setNames(c("0", "0.5", "1")), limits = 0:1) +  
  scale_y_continuous(breaks = (0:2)/2, limits = 0:1) +
  scale_color_discrete(name = "Condition\non") +
  scale_shape_discrete(name = "Condition\non")
ggsave("shareseq_condition_confounders.pdf", width = 6, height = 3.5)
ggsave("shareseq_condition_confounders.svg", width = 6, height = 3.5)


# Display KNN exchangeability results
knn_test %>%
  lapply(extract, c("prop_not_swapped", "p_value")) %>%
  sapply(unlist) %>%
  t %>%
  as.matrix %>%
  as.data.frame() %>%
  cbind(conditions) %>%
  subset(condition_on == "none" & error_mode == "none") %>%
  reshape2::melt(measure = c("prop_not_swapped", "p_value")) %>% 
  subset(variable %in% c("prop_not_swapped", "p_value")) %>%
  ggplot() +
  geom_point(aes(x = as.character(cell_count_cutoff), color = knockoff_type, y = value), stat = "identity", position = "dodge") +
  facet_wrap(~variable) +
  xlab("Cell count cutoff") + 
  geom_hline(aes(yintercept = value), data = data.frame(value = 0.5, variable = "prop_not_swapped")) + 
  ggtitle("KNN exchangeability diagnostic", "Real data")
ggsave(file = "KNN_exchangeability_test.pdf", width = 4, height = 2.5)
ggsave(file = "KNN_exchangeability_test.svg", width = 4, height = 2.5)

# Cell type composition
to_plot = set_up_share_skin_pseudobulk(conditions, 1)$pseudo_bulk_metadata %>%
  dplyr::mutate(
    ncell_exceeds_10  = total_cell_count >= 10,
    ncell_exceeds_100 = total_cell_count >= 100,
    ncell_exceeds_500 = total_cell_count >= 500
  ) %>%
  tidyr::pivot_longer(cols = dplyr::starts_with("ncell_exceeds"),
                      names_prefix = "ncell_exceeds_",
                      names_to = "cell_count_cutoff",
                      values_to = "passes_cutoff") %>%
  subset(passes_cutoff) %>%
  merge(supertypes, by.x = "largest_subpopulation", by.y = "celltype")
ggplot(to_plot) +
  geom_bar(aes(x = cell_count_cutoff, fill = supertype), position = "stack") +
  ylab("Number of metacells")

