library(pbivnorm)
predict_multivarite_probit_2 <- function(fit, formula, newdata, n_max = 1000) {
  model_matrix <- stats::model.matrix(formula, newdata)
  beta_dim <- fit@par_dims$beta
  omega_dim <- fit@par_dims$Omega
  stopifnot(NCOL(model_matrix) == beta_dim[2])
  stopifnot(omega_dim[1] == 2)
  beta_pars <- paste0("beta[", rep(1:beta_dim[1], each = beta_dim[2]), ",",
                      rep(1:beta_dim[2], beta_dim[1]), "]")
  Omega_pars <- paste0("Omega[",  rep(1:omega_dim[1], each = omega_dim[2]), ",",
                       rep(1:omega_dim[2], omega_dim[1]), "]")
  draws_df <- as.data.frame(rstan::extract(fit, pars = c(beta_pars, Omega_pars)))
  if (n_max < NROW(draws_df)) {
    i <- sample.int(NROW(draws_df), n_max, replace = FALSE)  
    draws_df <- draws_df[i, ]
  }
  
  res <- lapply (1:NROW(draws_df), function(j) {
    means <- sapply(1:beta_dim[1], function(i) {
      betas <- as.numeric(draws_df[j, paste0("beta.", i, ".", 1:beta_dim[2], ".")])
      model_matrix %*% betas  
    })
    corr <- draws_df[j, "Omega.1.2."]
    p11 <- pbivnorm(x = means[,1], y = means[,2], rho = corr)
    p10 <- pbivnorm(x = means[,1], y = -means[,2], rho = -corr)
    p01 <- pbivnorm(x = -means[,1], y = means[,2], rho = -corr)
    p00 <- 1 - (p11 + p10 + p01)
    p1_uncond <- p11 + p10
    p2_uncond <- p11 + p01
    p1_cond <- p11/(p11 + p01)
    p2_cond <- p11/(p11 + p10)# conditional prob of species 2 being present if species 1 is present
    p1_cond_0 <- p10/(p10 + p00)
    p2_cond_0 <- p01/(p01 + p00) # conditional prob of species 2 being present if species 1 is absent
    list(p1_uncond, p2_uncond, p1_cond, p2_cond, p1_cond_0, p2_cond_0)
  })
  list("p1_uncond" = sapply(res, "[[", 1),
       "p2_uncond" = sapply(res, "[[", 2),
       "p1_cond" = sapply(res, "[[", 3),
       "p2_cond" = sapply(res, "[[", 4),
       "p1_cond_0" = sapply(res, "[[", 5),
       "p2_cond_0" = sapply(res, "[[", 6))
}



summary(apply(pred_mul_un$p1_uncond, 1, mean))
summary(apply(pred_mul_un$p2_uncond, 1, mean))
summary(apply(pred_mul_un$p1_cond, 1, mean))
summary(apply(pred_mul_un$p2_cond, 1, mean))
summary(apply(pred_mul_un$p2_cond_0, 1, mean))
summary(apply(pred_mul_un$p1_cond_0, 1, mean))


