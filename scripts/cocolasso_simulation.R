# Test CoCoLASSO under a range of conditions ranging from ideal to realistic.
source("../scripts/setup.R")
conditions = read.table(
  header = T, 
  text = 
"knockoff_type condition_on cell_count_cutoff do_downsample
gaussian none  10 FALSE"
)
cat("\nReading in pseudo-bulk data...\n")
normalized_data = set_up_share_skin_pseudobulk(conditions, 1)

# For faster tests, use fewer TF's & targets
test_mode = T
if(test_mode){
  normalized_data$gene_metadata$is_tf = 
    rbinom(prob = 0.01, size = 1, n = length(normalized_data$gene_metadata$is_tf)) %>%
    as.logical()
  n_target = 2
} else {
  n_target = 50
}

# Each cluster is taken to be totally homogeneous. 
# Cells are measured with Gaussian error with the empirical cluster mean and some reasonable
# variance estimate -- maybe a trended common dispersion from edgeR. (For now, just sample variance.)
variance_true = normalized_data$pseudo_bulk_rna_var[which(normalized_data$gene_metadata$is_tf),] #in TPM^2
mean_true     = t(normalized_data$pseudo_bulk_rna)[ which(normalized_data$gene_metadata$is_tf),] #in TPM
pdf("cocolasso_simulation_snr.pdf")
plot(apply(mean_true, 1, var) %>% log10, xlab = "Variance across clusters",
     apply(variance_true, 1, mean) %>% log10, ylab = "Variance within clusters",
     main = "Signal to noise")
dev.off()
cluster_proportions = normalized_data$pseudo_bulk_metadata$total_cell_count %>% prop.table
simulate_one_cluster = function(cluster_idx, error_variance_scale = 1){
  n_cell = normalized_data$pseudo_bulk_metadata$total_cell_count[cluster_idx]
  n_tf = sum( normalized_data$gene_metadata$is_tf )
  simulated_cells =
    matrix( rnorm( n_tf*n_cell ), nrow = n_cell ) %*%
    Matrix::Diagonal(x = variance_true[,cluster_idx] %>% sqrt) 
    simulated_cells = simulated_cells*error_variance_scale
  
  simulated_cells %<>% sweep(2, mean_true[,cluster_idx], "+")
  simulated_cells
}
simulate_all_clusters = function(error_variance_scale){
  sim_expression_gaussian_mixture = matrix(0, 
                                           ncol = sum( normalized_data$gene_metadata$is_tf ), 
                                           nrow = sum(normalized_data$pseudo_bulk_metadata$total_cell_count) )
  position = 0
  for(cluster_idx in seq_along(normalized_data$pseudo_bulk_metadata$cluster)){
    cat(".")
    this_cluster = simulate_one_cluster(cluster_idx, error_variance_scale)
    sim_expression_gaussian_mixture[ position + (1:nrow(this_cluster)), ] = as.matrix(this_cluster)
    position = position + nrow(this_cluster)
  }
  sim_expression_gaussian_mixture
}
# In practice, measurement error dispersion is unknown. 
# We need to try CoCoLASSO with noisy estimates of it.
estimate_variance = function(sim_expression_gaussian_mixture){ 
  normalized_data$pseudo_bulk_metadata$position = cumsum( normalized_data$pseudo_bulk_metadata$total_cell_count )
  variance_estimates = sapply(
    seq_along(normalized_data$pseudo_bulk_metadata$cluster), 
    function(cluster_idx){
      cells_to_use = with(
        normalized_data$pseudo_bulk_metadata[cluster_idx,],
        (1:total_cell_count) + position - total_cell_count
      )
      apply(sim_expression_gaussian_mixture[cells_to_use,], 2, var)
    }
  )
  variance_estimates
}

# Run a single test
test_cocolasso = function(conditions, condition_idx){
  cat("\nTesting conditions:\n")
  write.table(conditions[condition_idx,], quote = F, sep = "\t")
  set.seed(conditions[condition_idx, "seed"])
  # Simulate TF's as perfectly homogeneous within each cluster
  cat("\nSimulating data...\n")
  sim_expression_gaussian_mixture = simulate_all_clusters(error_variance_scale = conditions[condition_idx, "error_variance_scale"] )
  cat("\nEstimating mean and variance... \n")
  # Estimate the per-gene, per-cluster variance
  variance_estimates = estimate_variance(sim_expression_gaussian_mixture)
  # Center & scale the features
  mean_x = mean_true %*% cluster_proportions
  var_x = sweep( mean_true, 1, mean_x, function(a, b) (a-b)^2) %*% cluster_proportions
  sim_expression_gaussian_mixture %<>% sweep(2,     mean_x,  "-")
  sim_expression_gaussian_mixture %<>% sweep(2, sqrt(var_x), "/")
  # Simulate targets
  cat("\nSimulating targets...\n")
  if(conditions[condition_idx, "response_type"]=="identity"){
    sim_targets_gaussian = simulate_all_clusters(error_variance_scale = conditions[condition_idx, "error_variance_scale"] )[, 1:n_target]
    true_active_sets = as.list(1:n_target)
  } else {
    logit = function(x) 1/(1 + exp(-4*x))
    FUN = switch(
      conditions[condition_idx, "response_type"], 
      linear=mean, 
      log_linear=prod, 
      threshold=function(x) all(x>0), 
      logit=function(x) sum(logit(x)),
      stop("response_type not recognized.\n")
    )
    sim_targets_gaussian = simulate_all_clusters(error_variance_scale = 0)[, 1:(2*n_target)]
    true_active_sets = list()
    for(target_idx in seq(n_target)){
      true_active_sets[[target_idx]] = 2*target_idx + c(-1:0)
      sim_targets_gaussian[,target_idx] = 
        sim_targets_gaussian[,true_active_sets[[target_idx]]] %>%
        apply(1, function(x) FUN(x)) %>%
        add(rnorm(nrow(sim_expression_gaussian_mixture)))
    } 
  }
  
  # Use either known or estimated measurement error variance, or pretend there is none.
  if(       conditions[condition_idx, "measurement_error"]=="true"){
    noise_covariance     =  conditions[condition_idx, "error_variance_scale"] * ( variance_true      %*% cluster_proportions ) / var_x 
  } else if(conditions[condition_idx, "measurement_error"]=="estimated"){
    noise_covariance     =                         ( variance_estimates %*% cluster_proportions ) / var_x 
  } else if(conditions[condition_idx, "measurement_error"]=="disregarded"){
    noise_covariance     =  0*(                     variance_estimates  %*% cluster_proportions )
  } else {
    stop("measurement_error mode not recognized.\n")
  }
  
  observed_gram_matrix = crossprod(sim_expression_gaussian_mixture) / nrow(sim_expression_gaussian_mixture)
  estimated_gram_matrix = observed_gram_matrix - Matrix::Diagonal( x = noise_covariance )
  closest_pd = BDcocolasso::ADMM_proj(as.matrix(estimated_gram_matrix))
  cocolasso_input_XX = closest_pd$mat*nrow(sim_expression_gaussian_mixture)
  # plot(diag(observed_gram_matrix), diag(estimated_gram_matrix), main = "Naive vs corrected feature variance")
  # test on many simulated targets
  cat("\nTesting COCOLASSO on each target...\n")
  early_precision = recall = precision = rep(NA, n_target)
  lambda = 10000000
  for(target_idx in 1:n_target){
    cat(".")
    cocolasso_input_XY = t(t(sim_targets_gaussian[,target_idx]) %*% sim_expression_gaussian_mixture)
    # Try lower lambda untin something comes up
    new_lambda = lambda
    iter = 0
    while(T){
      iter = iter + 1
      betahat = 
        BDcocolasso::lasso_covariance(
          p = ncol(sim_expression_gaussian_mixture), 
          n = nrow(sim_expression_gaussian_mixture),
          lambda = new_lambda, 
          penalty = "LASSO", 
          XX = cocolasso_input_XX, 
          Xy = cocolasso_input_XY, 
          beta.start = rep(0, ncol(observed_gram_matrix)) 
        )      
      if(iter>50 | max(abs(betahat$coefficients))>0){
        break
      }
      new_lambda = new_lambda / 2
    }
    dir.create(paste0("coefficients/condition_idx=", condition_idx), recursive = T)
    write.table(betahat$coefficients, paste0("coefficients/condition_idx=", condition_idx, "/target=", target_idx, ".csv"))
    n_correct     =    sum( which(betahat$coefficients!=0) %in% true_active_sets[[target_idx]] )
    n_discoveries = length( which(betahat$coefficients!=0) )
    n_true        = length( true_active_sets[[target_idx]] )
    
    early_precision[target_idx] = which.max(betahat$coefficients) %in% true_active_sets[[target_idx]]
    recall[target_idx] = n_correct / n_true
    precision[target_idx] = n_correct / n_discoveries
  }
  return( 
    data.frame( 
      early_precision = mean(early_precision),
      recall          = mean(recall),
      precision       = mean(precision) 
    ) 
  )
}


# We test several different settings.
conditions = read.table(
  header = T, 
  text = 
    "response_type measurement_error seed condition_on error_variance_scale
          identity              true    1         none                  0.3
          identity       disregarded    1         none                  0.3
          identity              true    1         none                  0.2
          identity       disregarded    1         none                  0.2
          identity              true    1         none                  0.1
          identity       disregarded    1         none                  0.1
          identity              true    1         none                  0.05
          identity       disregarded    1         none                  0.05
          identity         estimated    1         none                  0.5
          identity       disregarded    1         none                  0.5
            linear         estimated    1         none                  0.5
        log_linear         estimated    1         none                  0.5
         threshold         estimated    1         none                  0.5
             logit         estimated    1         none                  0.5
  "
)[1:8,]
results = lapply(
  seq_along(conditions[[1]]), 
  function(condition_idx){
    test_cocolasso(conditions, condition_idx)
  }
) 
print(results)
saveRDS(results, "cocolasso_simulation.Robj")
results %>%
  data.table::rbindlist() %>% 
  cbind(conditions) %>%
  write.csv("cocolasso_simulation.csv")

results %>%
  data.table::rbindlist() %>% 
  cbind(conditions) %>% 
  write.csv()
