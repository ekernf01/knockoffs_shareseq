source("../scripts/setup.R")
sanitize_names = function(df){
  colnames(df) = colnames(df) %>% 
    gsub("empirical_fdr", "observed_fdr", .)  %>% 
    gsub("nominal_fdr", "expected_fdr", .) 
  return(df)
}
ggplot2::theme_update(text = element_text(family = "ArialMT"))
conditions = read.csv("experiments_to_run.csv", row.names = 1)
# Plot big calibration complicated many facets such oof wow
all_calibration = knn_test = list()
for(condition_idx in seq(nrow(conditions))){
  # Set up condition-specific environment
  attach(conditions[condition_idx,], warn.conflicts = F)
  cat("\nRetrieving condition:\n")
  write.table(t(conditions[condition_idx,]), quote = F, sep = "\t")
  new_wd = create_experiment_path(conditions, condition_idx)
  normalized_data = set_up_share_skin_pseudobulk(conditions, condition_idx)
  # These help condense results from rlookc::simulateY()
  summarize_simulation = function(X) {
    X %>%
      colMeans %>%
      (function(x) data.frame(expected_fdr = gsub("^X", "", names(x)) %>% as.numeric, observed_fdr = x)) %>%
      dplyr::mutate( num_discoveries  = NA,  moe_95 = NA)
  }
  read_simulation_if_exists = function(error_mode){
    if(file.exists("calibration/average_case_calibration.csv")){
      read.csv("calibration/average_case_calibration.csv", row.names = 1) %>%
        summarize_simulation %>%
        cbind(evaluation=paste0("simulated y with error_mode=", error_mode)) %>%
        return
    } else {
      return(NULL)
    }
  }
  # Run KNN exchangeability diagnostic
  # Set n_neighbors low because there's not many observations
  withr::with_dir(new_wd, {
    knn_test[[paste0(condition_idx, "_")]] = NULL
    try( silent = T, {
      knn_test[[paste0(condition_idx, "_")]] = rlookc::KNNTest(X = readRDS("output_knockoffs/original_features.Rda"),
                                                X_k = readRDS("output_knockoffs/knockoffs.Rda"),
                                                n_neighbors = 5)
    })
  })
  # Grab outputs and make them uniform
  withr::with_dir(new_wd, {
    all_calibration[[condition_idx]] = list(
      read.csv("calibration/decoy_calibration.csv", row.names = 1) %>%
        cbind(evaluation="decoy") %>% 
        sanitize_names,
      read.csv("calibration/chip_calibration.csv", row.names = 1) %>%
        cbind(evaluation="chip") %>% 
        sanitize_names,
      read_simulation_if_exists(error_mode=error_mode)
    )  %>% 
      (function(X) (X[!sapply(X, is.null)])) %>%
      lapply(extract, c("expected_fdr", "observed_fdr", "evaluation")) %>%
      Reduce(f=rbind) %>%
      merge(conditions[condition_idx,])
  })
}
all_calibration %<>% Reduce(f = rbind)
all_calibration$celltype %<>% factor(levels = c("skin", "keratinocyte", "pbmc", "tcell", "pbmc_subset"))
all_calibration$tf_activity_type %<>% factor(levels = c("rna", "motif", "both"))

# These plots are for an internal talk
create_experiment_path(conditions, 49+10) %>% 
  file.path("all_hypotheses.csv.gz") %>%
  read.csv() %>% 
  subset(q <= 0.1, select = "Gene1") %>% 
  table %>%
  table %>%
  as.data.frame() %>%
  set_colnames(c("outdegree_char", "frequency")) %>%
  dplyr::mutate(outdegree = as.numeric(outdegree_char)) %>%
  ggplot() + geom_point(
    aes(y = frequency, x = outdegree)
  ) + 
  scale_x_log10() + 
  scale_y_log10() +
  xlab("out-degree") +
  ylab("frequency") + 
  ggtitle("All PBMC")
create_experiment_path(conditions, 65+10) %>% 
  file.path("all_hypotheses.csv.gz") %>%
  read.csv() %>% 
  subset(q <= 0.1, select = "Gene1") %>% 
  table %>%
  table %>%
  as.data.frame() %>%
  set_colnames(c("outdegree_char", "frequency")) %>%
  dplyr::mutate(outdegree = as.numeric(outdegree_char)) %>%
  ggplot() + geom_point(
    aes(y = frequency, x = outdegree)
  ) + 
  scale_x_log10() + 
  scale_y_log10() +
  xlab("out-degree") +
  ylab("frequency") + 
  ggtitle("T cell")


# Extract a few key numbers mentioned in the text (number of discoveries in certain analyses)
total_findings_1_percent = list(
  skin_basic = 
    create_experiment_path(conditions, 1+10) %>% 
    file.path("all_hypotheses.csv.gz") %>%
    read.csv() %>% 
    subset(q <= 0.1) %>% 
    nrow,
  pbmc_basic = 
    create_experiment_path(conditions, 49+10) %>% 
    file.path("all_hypotheses.csv.gz") %>%
    read.csv() %>% 
    subset(q <= 0.1) %>% 
    nrow, 
  skin_100 = 
    create_experiment_path(conditions, 2+10) %>% 
    file.path("all_hypotheses.csv.gz") %>%
    read.csv() %>% 
    subset(q <= 0.1) %>% 
    nrow,
  pbmc_100 = 
    create_experiment_path(conditions, 50+10) %>% 
    file.path("all_hypotheses.csv.gz") %>%
    read.csv() %>% 
    subset(q <= 0.1) %>% 
    nrow
)
total_findings_1_percent

# Print the motif network summary stats
for(i in c(73, 32)){
  x = create_experiment_path(conditions, i) %>% 
    file.path("all_hypotheses.csv.gz") %>%
    read.csv()
  subset(x, is_testable)$Gene1 %>% unique %>% length %>% print
  subset(x, is_testable)$Gene2 %>% unique %>% length %>% print
  subset(x, is_testable)[c("Gene1", "Gene2")] %>% apply(1, paste, collapse = " ") %>% unique %>% length %>% print
}

# Plot naive versus "shrinkage" (corpcor gaussian) knockoffs
all_calibration %>%
  subset(
    T & 
      knockoff_type %in% c("gaussian", "permuted") &
      celltype %in% c("skin", "pbmc") &
      tf_activity_type == "rna"  &
      condition_on == "none"     & 
      error_mode == "none"       &            
      seed==1                    &
      require_motif_support == F & 
      only_motif_support == F    &
      evaluation == "simulated y with error_mode=none"
    ) %>%
  ggplot() +
  geom_abline(aes(slope = 1, intercept = 0)) +
  geom_point(aes(x = expected_fdr, y = observed_fdr, color = knockoff_type, shape = knockoff_type)) +
  facet_grid(celltype ~ cell_count_cutoff) +
  ggtitle("Conditional independence testing on multi-omics data", "Real TF expression and simulated targets") + 
  theme(text = element_text(family = "ArialMT")) +  
  scale_x_continuous(breaks = (0:2)/2, limits = 0:1) +  
  scale_y_continuous(breaks = (0:2)/2, limits = 0:1) + 
  xlab("Reported FDR") + 
  ylab("True FDR")
ggsave("shareseq_naive_vs_gaussian.pdf", width = 5, height = 4)
ggsave("shareseq_naive_vs_gaussian.svg", width = 5, height = 4)

# A supplemental table
all_calibration %>%
  subset(expected_fdr==0.2) %>%
  write.csv("results_fdr=20pct.csv")

# Plot initial ChIP results
all_calibration %>%
  subset(
    T & 
      knockoff_type %in% c("gaussian") &
      celltype %in% c("skin", "keratinocyte", "pbmc", "tcell") &
      tf_activity_type == "rna"  &
      condition_on == "none"     & 
      # cell_count_cutoff == 10    & 
      error_mode == "none"       &            
      seed==1                    &
      require_motif_support == F & 
      only_motif_support == F    &
      grepl("chip", evaluation)  
  ) %>% 
  ggplot() +
  geom_abline(aes(slope = 1, intercept = 0)) +
  geom_point(aes(
    x = expected_fdr, 
    y = observed_fdr,
    color = as.character(cell_count_cutoff) )) +
  ggtitle("TRN inference on multi-omics data", "Real data") + 
  theme(text = element_text(family = "ArialMT")) +  
  scale_x_continuous(breaks = (0:2)/2, limits = 0:1) +  
  scale_y_continuous(breaks = (0:2)/2, limits = 0:1) + 
  facet_grid(~celltype) + 
  labs(color = "Cell count cutoff") + 
  xlab("Reported FDR") + 
  ylab("True FDR") + 
  theme(legend.position = "bottom")
ggsave("shareseq_cellcount_cutoff.pdf", width = 4, height = 2.5)
ggsave("shareseq_cellcount_cutoff.svg", width = 4, height = 2.5)

# Plot initial decoy results
all_calibration %>%
  subset(
    T & 
      knockoff_type %in% c("gaussian") &
      celltype %in% c("skin", "keratinocyte", "pbmc", "tcell") &
      tf_activity_type == "rna"  &
      condition_on == "none"     & 
      # cell_count_cutoff == 10    & 
      error_mode == "none"       &            
      seed==1                    &
      require_motif_support == F & 
      only_motif_support == F    &
      grepl("decoy", evaluation)  
  ) %>% 
  ggplot() +
  geom_abline(aes(slope = 1, intercept = 0)) +
  geom_point(aes(
    x = expected_fdr, 
    y = observed_fdr,
    color = as.character(cell_count_cutoff) )) +
  ggtitle("TRN inference on multi-omics data", "Real data") + 
  theme(text = element_text(family = "ArialMT")) +  
  scale_x_continuous(breaks = (0:2)/2, limits = 0:1) +  
  scale_y_continuous(breaks = (0:2)/2, limits = 0:1) + 
  facet_grid(~celltype) + 
  labs(color = "Cell count cutoff") + 
  xlab("Reported FDR") + 
  ylab("Observed proportion of decoys") + 
  theme(legend.position = "bottom") +
  geom_hline(aes(yintercept = 100/1078))
ggsave("shareseq_cellcount_cutoff_decoys.pdf", width = 4, height = 2.5)
ggsave("shareseq_cellcount_cutoff_decoys.svg", width = 4, height = 2.5)

# Plot results with/without error in X
all_calibration %>%
  subset(
    T & 
      knockoff_type %in% c("gaussian") &
      celltype %in% c("skin", "keratinocyte", "pbmc", "tcell") &
      tf_activity_type == "rna"  &
      condition_on == "none"     & 
      error_mode == "resample" &
      # cell_count_cutoff == 500    &
      # error_mode == "none"       &            
      seed==1                    &
      require_motif_support == F & 
      only_motif_support == F    &
      grepl("simulated", evaluation)
  ) %>% 
  ggplot() +
  geom_abline(aes(slope = 1, intercept = 0)) +
  geom_point(aes(
    x = expected_fdr, 
    y = observed_fdr,
    color = as.character(cell_count_cutoff) )) +
  facet_grid(1~celltype) +
  ggtitle("Conditional independence testing on multi-omics data", "Real TF expression with additional error and simulated targets") + 
  theme(text = element_text(family = "ArialMT")) +  
  scale_x_continuous(breaks = (0:2)/2, limits = 0:1) +  
  scale_y_continuous(breaks = (0:2)/2, limits = 0:1) + 
  labs(color = "Error mode") + 
  xlab("Reported FDR") + 
  ylab("True FDR") + 
  labs(color = "Cell count cutoff") + 
  theme(legend.position = "bottom")
ggsave("shareseq_measurement_error.pdf", width = 5.5, height = 3)
ggsave("shareseq_measurement_error.svg", width = 5.5, height = 3)

# Show results with different TF activities
all_calibration %>%
  subset(
    T & 
      knockoff_type %in% c("gaussian") &
      celltype %in% c("skin", "keratinocyte", "pbmc", "tcell") &
      condition_on == "none"     & 
      cell_count_cutoff == 10    &
      # celltype %in% c("pbmc") &
      error_mode == "none"       &
      seed==1                    &
      require_motif_support == F & 
      only_motif_support == F    &
      grepl("chip", evaluation)  == T
  ) %>% 
  ggplot() +
  geom_abline(aes(slope = 1, intercept = 0)) +
  geom_point(aes(x = expected_fdr, y = observed_fdr, color = tf_activity_type)) +
  facet_grid( ~ celltype ) +
  ggtitle("Calibration by type of TF activity", "Real data") + 
  theme(text = element_text(family = "ArialMT")) +  
  scale_x_continuous(breaks = ((0:2)/2) %>% setNames(c("0", "0.5", "1")), limits = 0:1) +  
  scale_y_continuous(breaks = (0:2)/2, limits = 0:1)  + 
  xlab("Reported FDR") + 
  ylab("True FDR") + 
  theme(legend.position = "bottom") + 
  scale_color_manual(values = c("red", "blue", "purple"))
ggsave("shareseq_motif_scores.svg", width = 4.5, height = 2.5)
ggsave("shareseq_motif_scores.pdf", width = 4.5, height = 2.5)

# Show results with/without conditioning on confounders
all_calibration %>%
  subset(
    T & 
      knockoff_type %in% c("gaussian") &
      celltype %in% c("skin", "keratinocyte", "pbmc", "tcell") &
      tf_activity_type == "both"  &
      error_mode == "none"       &
      seed==1                    &
      cell_count_cutoff == 10   & 
      require_motif_support == F &
      only_motif_support == F    &
      grepl("chip", evaluation)  == T
  ) %>% 
  ggplot() +
  geom_abline(aes(slope = 1, intercept = 0)) +
  geom_point(aes(
    x = expected_fdr, 
    y = observed_fdr,
    color = condition_on )) +
  facet_grid(1~celltype) +
  ggtitle("TRN inference on multi-omics data", "Real data") + 
  theme(text = element_text(family = "ArialMT")) +  
  scale_x_continuous(breaks = (0:2)/2, limits = 0:1) +  
  scale_y_continuous(breaks = (0:2)/2, limits = 0:1) + 
  xlab("Reported FDR") + 
  ylab("True FDR") + 
  theme(legend.position = "bottom") + 
  scale_color_manual(values = c("green", "black"))
ggsave("shareseq_condition_confounders.pdf", width = 5.5, height = 3)
ggsave("shareseq_condition_confounders.svg", width = 5.5, height = 3)

# Show results with motif support
all_calibration %>%
  subset(
    T & 
      condition_on == "pca" &
      celltype %in% c("skin", "keratinocyte", "pbmc", "tcell") &
      cell_count_cutoff == 10 &
      tf_activity_type=="both" &
      grepl("chip", evaluation)  == T 
  ) %>% 
  dplyr::mutate(
    hypothesis_screening = 
      ifelse(
        only_motif_support, 
        "motif only", 
        ifelse(
          require_motif_support,
          "motif and knockoff-based",
          "knockoff-based"
        )
      ) 
    ) %>%
  ggplot() +
  geom_abline(aes(slope = 1, intercept = 0)) +
  geom_point(aes(
    x = expected_fdr, 
    y = observed_fdr,
    color = hypothesis_screening )) +
  labs(color = "Analysis method") +
  facet_grid(~celltype) +
  ggtitle("TRN inference on multi-omics data", "Real data") + 
  theme(text = element_text(family = "ArialMT")) +  
  scale_x_continuous(breaks = (0:2)/2, limits = 0:1) +  
  scale_y_continuous(breaks = (0:2)/2, limits = 0:1) + 
  xlab("Reported FDR") + 
  ylab("True FDR") + 
  theme(legend.position = "bottom") + 
  scale_color_manual(values = c("red", "orange", "yellow"))
ggsave("shareseq_motif_support.pdf", width = 6, height = 3.5)
ggsave("shareseq_motif_support.svg", width = 6, height = 3.5)

# Display KNN exchangeability results
knn_test %>%
  lapply(extract, c("prop_not_swapped", "p_value")) %>%
  sapply(unlist) %>%
  t %>%
  as.matrix %>%
  as.data.frame() %>%
  tibble::rownames_to_column("condition") %>%
  merge(
    conditions %>% tibble::rownames_to_column("condition") %>%dplyr::mutate(condition = paste0(condition, "_")) ,
    by = "condition"
  ) %>%
  subset(condition_on == "none" & error_mode == "none" & tf_activity_type == "rna" & celltype %in% c("skin", "pbmc")) %>%
  reshape2::melt(measure = c("prop_not_swapped", "p_value")) %>% 
  subset(variable %in% c("prop_not_swapped", "p_value")) %>%
  ggplot() +
  geom_point(aes(x = as.character(cell_count_cutoff), color = knockoff_type, y = value), stat = "identity", position = "dodge") +
  facet_grid(variable ~ celltype) +
  xlab("Cell count cutoff") + 
  geom_hline(aes(yintercept = value), data = data.frame(value = 0.5, variable = "prop_not_swapped")) + 
  ggtitle("KNN exchangeability diagnostic", "Real data")
ggsave(file = "KNN_exchangeability_test.pdf", width = 4, height = 4)
ggsave(file = "KNN_exchangeability_test.svg", width = 4, height = 4)

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

