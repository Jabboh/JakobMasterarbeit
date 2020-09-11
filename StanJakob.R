##################################Analysis of Culex perexiguus & Anopheles troparvus:
#Comparing Univariate probit models with multivariate probit models. Procedure:
#1. Data preparation
#2. Fitting the most complex models (all environmental covariates 
#+ according interaction and quadratic terms): a. probit Culex, b. Probit Anopheles
#c. multivariate probit Anopheles & Culex
#3. Internal Model validation of these models by analyzing its residuals with DHARMa
#4. Selecting the "best" model by using WAIC on the univariate probit models
#>> Finding the "most appropriate" specifications of covariates
#5. Redoing the internal validation for the "most appropriate" model.
#6. Results
#6.1. In-sample (not on test data) or maybe model exploration?:
#a. Comparing the coefficients: Size and credibility intervals
#b. Correlations between the responses / Residual Correlation Parameter of gjam
#c. Response Curves 
#d. Variable importance
#6.2. Out-of-Sample
#a. Conditional Predictions of gjam vs. Unconditional Predictions of 
#gjam vs. predictions of univariate model
#b. Comparing the uncertainty of the different "prediction"-types
rm(list=ls())
library(here)
#install.packages("")

#loading packages
library(readxl) #to read in the data
library(rstan)
library(rstanarm) #Doing Bayesian probit GLMs for single species
library(mvtnorm)
library(pROC) #calculating the AUC
# library(ggplot2) #for plotting
# library(DHARMa) # for checking in-sample validity of the models
library(dplyr) # for simplified syntax and neater code
library(bayesplot) #Some handy features for plotting in the Bayesian realm
# library(psych)
#### 1.Data Preparation
#read in the data (Monthly species PA data for seven different mosquito species
#and according environmental covariates)
df <- read_excel(here("Data/MonthlyData.xlsx"))

#Checking the rough structure of the data set
# summary(df$Fecha) # Dates in 2006 dont make any sense. I assume that they put by
#accident 2006 instead of 2010. So I change these dates
df$Fecha <- as.POSIXct(sub("2006", "2010", df$Fecha))

#Transform "Mes" (month when record was taken) into a factor variable
df$Mes <- factor(df$Mes, levels =c("Abril", "Mayo", "Junio", "Julio", "Agosto",
                                   "Septiembre"))

#Reading in the spatial coordinates of the different trap locations
coords <- read_excel(here("Data/Traps_coordenadas_geograficas.xls"))

#Trap and Area mean the same thing, so we change the name in df from "area" to "trap"
names(df)[names(df)=="Area"] <- "trap"
#Canada is spelt differently, so I change the spelling in df to "Cañada"
df[,"trap"] <- lapply(df[,"trap"], gsub, pattern = "Ca?da", replacement = "Cañada",
                      fixed = T)
#adding lon-lat column to the data frame df
df <- merge(df, coords[, c("trap", "Norte", "Oeste")], by="trap", all.x= T, sort = F)

#Selecting our two species (perexiguus & Anopheles troparvus) for the analysis. Our
#rationale is to select the two species that roughly occur in 50 % of the observations
#in order to get the maximum variation. We ignore the column "An_atroparvus", 
#because it does not seem to be PA-data (rather abundance data).

#extract the PA data for all species
spec <- df[, 7:14]
#deleting the An_atroparvus column
spec[,"An_atroparvus"] <- NULL
#Taking a look at occurence rate across all observations (the Mean)
# summary(spec)
#Hence, we select Cxperpre and Anatropre for our first analysis:
y <- as_tibble(spec[,c("Cxperpre", "Anatropre")])

#Normalizing Covariates: Otherwise, interaction terms hardy interpretable and skewed
df[,17:36] <- scale(df[,17:36])

#Split the data set in training (70 %) and test (30%) set
train_size <- floor(0.7 * nrow(y))
#set a random seed for replicability
set.seed(333)
#sample the training IDs
train_id <- sample(seq_len(nrow(y)), size = train_size)

#partition data into train and test set
train <- df[train_id, ]
test <- df[-train_id, ]
y_train <- y[train_id, ]
y_test <- y[-train_id, ]

#' Use multiple cores in stan, and only recompile models when stan code has changed
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#' Create data sets for stan multivariate probit model
make_multivariate_data <- function(formula, ydata, data) {
  model_matrix <- stats::model.matrix(formula, data)
  list(
    K = NCOL(model_matrix),
    D = NCOL(ydata),
    N = NROW(data),
    y = as.matrix(ydata),
    x = model_matrix
  )
}


#' Model formula - only the right hand side (i.e. the side for the predictors)
rhs_formula <- ~ Mes + IA_500 + NDVIBEF_2000 + I(IA_500^2) + I(NDVIBEF_2000^2)

#' Create data sets for the different stan models
ap_cp_data <- make_multivariate_data(rhs_formula, ydata = y_train, data = train)
ap_data <- make_multivariate_data(rhs_formula, ydata = y_train[,"Anatropre", drop = FALSE], 
                                  data = train)
cp_data <- make_multivariate_data(rhs_formula, ydata = y_train[,"Cxperpre", drop = FALSE], 
                                  data = train)


#' Fit the CP with multiivariate probit model (Wenn die Ergebnisse Sinn ergeben, werde ich Rstan
#' Version benutzen (weil ich hab ja DHARMa und so alles mit rstan gemacht...))

cp_fit <- stan(
  file = "multivariate_probit.stan",
  data = cp_data, 
  chains = 4, 
  iter = 4000,
  refresh = 0,
  seed = 2
)

#THe Rhat-statistic measures the ratio of the average variance of draws within each chain to 
#the variance of the pooled draws across chains; if all chains are at equilibrium, these will 
#be the same and R^ will be one. If the chains have not converged to a common distribution, the
#R^ statistic will be greater than one (see Gelman et al. 2013, Stan Development Team 2018).
#The rule of thumb is that R-hat values for all less than 1.1 

maxRhat <- function(fit) {
  max(summary(fit)$summary[, "Rhat"], na.rm = TRUE)  
}

maxRhat(cp_fit)

system.time(ap_fit <- stan(
  file = "multivariate_probit.stan",  # Stan program
  data = ap_data,    # named list of data
  chains = 4,             # number of Markov chains
  iter = 4000,            # total number of iterations per chain
  refresh = 0,             # no progress shown
  seed = 2
))

maxRhat(ap_fit) #R-hat compares

ap_cp_fit <- stan(
  file = "multivariate_probit.stan",  # Stan program
  data = ap_cp_data,    # named list of data
  chains = 4,             # number of Markov chains
  iter = 4000,            # total number of iterations per chain
  refresh = 0,             # no progress shown
  seed = 2
)

maxRhat(ap_cp_fit)


extract_betas <- function(fit, species = 1) {
  beta_dim <- fit@par_dims$beta
  stopifnot(species <= beta_dim[1])
  beta_pars <- paste0("beta[", species, ",", 1:beta_dim[2], "]")
  as.data.frame(rstan::extract(fit, pars = beta_pars))
}

#' # Comparing betas

apcp_cp <- extract_betas(ap_cp_fit, 1)
cp_cp <- extract_betas(cp_fit, 1)



plot_betas <- function(x, y, x_label, y_label, title_species = "") {
  title <- paste(x_label, "vs.", y_label, "Coefficients and their 95 % - Credibility Intervals \n for", title_species)
  posterior_1 <- as.matrix(x) #rstanarm nimmt den median als Schätzer für den Parameter
  posterior_2 <- as.matrix(y)
  # colnames(posterior_2) <- colnames(posterior_1) <- colnames(Anatropre_Cxperpre_data$x)
  combined <- rbind(mcmc_intervals_data(posterior_1, prob_outer = .95), 
                    mcmc_intervals_data(posterior_2, prob_outer = .95))
  combined$model <- rep(c(x_label, y_label), each = ncol(posterior_1))
  #prob_outer defines the credibility interval in our plot (here .975)
  theme_set(bayesplot::theme_default())
  pos <- position_nudge(y = ifelse(combined$model == x_label, 0, 0.1))
  coef_cp <- ggplot(combined, aes(x = m, y = parameter, color = model)) + 
    geom_point(position = pos) +
    geom_vline(xintercept = 0, linetype="dotted", color = "black", size=.5) +
    geom_errorbar(aes(xmin = ll, xmax = hh), position = pos, width = .1) +
    ggtitle(title) + 
    xlab("Value") + ylab("Coefficient") + labs(color="Model") 
  coef_cp
}

colnames(cp_cp) <- colnames(apcp_cp) <- colnames(cp_data$x)
plot_betas(cp_cp, apcp_cp, "Univariate", "Multivariate", "Culex perexiguus")

#' Correlations between posterior coefficient estimates
plot(cor(cp_cp), cor(apcp_cp))
abline(0, 1)

#' # Comparison of predictions
#' 

# The following code can probably be made much faster using the pbivnorm from the package of the same name
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
    corr <- diag(omega_dim[1])
    corr[lower.tri(corr)] <- corr[upper.tri(corr)] <- draws_df[j, "Omega.1.2."]
    p11 <- apply(means, 1, function(m) mvtnorm::pmvnorm(lower=0, upper = Inf, mean = m, corr = corr))
    p10 <- apply(means, 1, function(m) mvtnorm::pmvnorm(lower=c(0, -Inf), upper = c(Inf, 0), mean = m, corr = corr))
    p01 <- apply(means, 1, function(m) mvtnorm::pmvnorm(lower=c(-Inf, 0), upper = c(0, Inf), mean = m, corr = corr))
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

set.seed(2)
system.time(ap_cp_multi_p <- predict_multivarite_probit_2(ap_cp_fit, rhs_formula, train))

#Bin Hier!!!!!!!!

predict_multivarite_probit_1 <- function(fit, formula, newdata, n_max = 4000, type = "link") {
  model_matrix <- stats::model.matrix(formula, newdata)
  beta_dim <- fit@par_dims$beta
  omega_dim <- fit@par_dims$Omega
  stopifnot(NCOL(model_matrix) == beta_dim[2])
  stopifnot(omega_dim[1] == 1, beta_dim[1] == 1)
  beta_pars <- paste0("beta[1,", 1:beta_dim[2], "]")
  draws_df <- as.data.frame(rstan::extract(fit, pars = beta_pars))
  if (n_max < NROW(draws_df)) {
    i <- sample.int(NROW(draws_df), n_max, replace = FALSE)  
    draws_df <- draws_df[i, ]
  }
  means <- apply(draws_df, 1, function(betas) model_matrix %*% betas)
  if (type == "link") return(means)
  if (type == "response") return(pnorm(means)) #Ich glaube hier ist was fishy!
}

#also type = response funktioniert meiner Meinung nicht
cp_p <- predict_multivarite_probit_1(cp_fit, rhs_formula, train, type = "response")



cp_multi_p <- ap_cp_multi_p$p1_uncond

#' rstanarm version of probit regression
full_formula <- as.formula(paste("Cxperpre ~", as.character(rhs_formula)[2]))
cp_arm_fit <- stan_glm(full_formula, data = train, refresh = 0, family = binomial(link = "probit"), 
                       init_r = 1.5, seed = 333)
cp_arm_p <- t(posterior_linpred(cp_arm_fit, seed = 23, draws = 4000, transform = TRUE))

#' Mean posterior predictions are the same for the 325 plots
op <- par(no.readonly = TRUE)


plot(rowMeans(cp_p), rowMeans(cp_multi_p), asp = 1, xlim = c(0, 1), ylim = c(0, 1))
abline(0, 1)

#with rstanarm
plot(rowMeans(cp_arm_p), rowMeans(cp_multi_p), asp = 1, xlim = c(0, 1), ylim = c(0, 1))
abline(0, 1)

par(op)
#multivariate and univariate models are more or less identical



#' # Uncertainties for the predicted plots - here expressed as the sd of posterior probability of occurrence
op <- par(no.readonly = TRUE)
par(mfrow = c(2,2))
abline(0, 1)
plot(apply(cp_p, 1, sd), apply(cp_multi_p, 1, sd), asp = 1)
abline(0, 1)
par(op)

#with rstanarm
plot(apply(cp_arm_p, 1, sd), apply(cp_multi_p, 1, sd), asp = 1)
abline(0, 1)


#multivariate and univariate models are more or less identical
############################BIN HIER
#' Conditional predictions
cp_multi_cond_p <- ap_cp_multi_p$p1_cond
cp_multi_cond_0_p <- ap_cp_multi_p$p1_cond_0

op <- par(no.readonly = TRUE)
par(mfrow = c(1, 2))
plot(rowMeans(cp_p), rowMeans(cp_multi_cond_p), col = "red", xlim = c(0, 1), ylim = c(0, 1),
     xlab = "Univariate prediction", ylab = "Conditional prediction",
     main = "All sites predicted twice\nother species assumed present or absent")
points(rowMeans(cp_p), rowMeans(cp_multi_cond_0_p), col = "blue")
abline(0, 1)

ap_1 <- train$Anatropre == 1
ap_0 <- train$Anatropre == 0

plot(rowMeans(cp_p)[ap_1], rowMeans(cp_multi_cond_p)[ap_1], col = "red", xlim = c(0, 1), ylim = c(0, 1),
     xlab = "Univariate prediction", ylab = "Conditional prediction",
     main = "Conditioned on realised \npresence of other species")
points(rowMeans(cp_p)[ap_0], rowMeans(cp_multi_cond_0_p)[ap_0], col = "blue")
abline(0, 1)
par(op)

#' AUC of predictions: conditional vs. unconditional
auc_uncond <- roc(train$Cxperpre, rowMeans(cp_multi_p))
auc_uncond$auc # 0.8
auc_cond <- roc(train$Cxperpre, ifelse(train$Anatropre == 1, rowMeans(cp_multi_cond_p), rowMeans(cp_multi_cond_0_p)))
auc_cond$auc # 0.817
