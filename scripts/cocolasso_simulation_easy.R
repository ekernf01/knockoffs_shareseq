compound_symmetric = matrix(0.5, nrow = 100, ncol = 100) + diag(x = rep(0.5, 100))
R = chol(compound_symmetric)
do_one = function(noise_scale){
  X     = rnorm(1e5) %>% matrix(ncol = 100) %>% multiply_by_matrix(R)
  Z = X + rnorm(1e5) %>% matrix(ncol = 100)*noise_scale
  Y = rowMeans(X[, 1:5] - X[,6:10])  + rnorm(1e3) %>% as.matrix() 
  noise_variance = diag(rep(1, ncol(Z)))*(noise_scale^2)
  beta_cocolasso = BDcocolasso::lasso_covariance(
    Xy = t(Z) %*% Y, XX = BDcocolasso::ADMM_proj(cov(Z) - noise_variance)$mat,
    n = 100, p = 100, penalty = "LASSO", beta.start = rep(0, ncol(Z)), lambda = 200
  )
  beta_lasso = BDcocolasso::lasso_covariance(
    Xy = t(Z) %*% Y, XX = cov(Z),
    n = 100, p = 100, penalty = "LASSO", beta.start = rep(0, ncol(Z)), lambda = 200
  )
  return(data.frame(
    true = c(rep(1, 5), rep(-1, 5), rep(0, 90)) %>% as.character(), 
    rank_lasso     = rank( abs( beta_lasso$coefficients  )), 
    rank_cocolasso = rank( abs( beta_cocolasso$coefficients )), 
    lasso     =  beta_lasso$coefficients , 
    cocolasso =  beta_cocolasso$coefficients , 
    noise_scale = noise_scale
  ))
}
df = 
  lapply(rep(3, 10), do_one) %>%
  data.table::rbindlist()
ggplot(subset(df)) + 
  geom_point(aes((cocolasso), (lasso), colour = true)) + 
  facet_wrap(~noise_scale) + 
  ggtitle("CoComparison")

with(df, table(detected_by_lasso     = rank_lasso>90, 
               detected_by_cocolasso = rank_cocolasso>90, 
               truth = true)) %>%
  View

# Lok at the Gram matrix
plot(cov(X), cov(Z), col = "red", pch = "*", 
     xlab ="true", ylab = "estimate", main = "Gram matrix entries")
points(cov(X), cov(Z) - noise_variance, col = "blue", pch = "-")
abline(a=0, b = 1)
legend(x = 0.4,
       y = 9, 
       col = c("red", "blue", "black"), 
       legend = c("naive", "unbiased", "y=x"), 
       pch = c("*", "-", "l"))
