###############################Results of my analyis#############################################
#1. Data preparation and fitting of final model
#6. Results
#6.1. In-sample (not on test data):
#a. Comparing the coefficients: Size and credibility intervals
#b. Correlations between the responses / Residual Correlation Parameter of multivariate model
#c. Response Curves 
#d. Variable importance
#6.2. Out-of-Sample
#a. Conditional Predictions of multivariate vs. Unconditional Predictions of 
#univariate vs. predictions of multivariate model
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
library(pbivnorm) #To calculate CDF of bivariate normal distribution
library(DALEX)
library(matrixStats) #to calculate quantiles
##########################################1.Data Preperation#######################################
#read in the data (Monthly species PA data for seven different mosquito species
#and according environmental covariates)
df <- read_excel(here("Data/MonthlyData.xlsx"))

#Checking the rough structure of the data set
# summary(df$Fecha) # Dates in 2006 dont make any sense. I assume that they put by
#accident 2006 instead of 2010. So I change these dates
df$Fecha <- as.POSIXct(sub("2006", "2010", df$Fecha))

#Transform "Month" (month when record was taken) into a factor variable
df$Mes <- factor(df$Mes, levels =c("Abril", "Mayo", "Junio", "Julio", "Agosto",
                                   "Septiembre"))
#Rename months in English
names(df)[names(df) == "Mes"] <- "Month"
levels(df$Month) <- c("April","May","June", "July", "August", "September")

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
rhs_formula <- ~ Month + IA_500 + NDVIBEF_2000 + I(IA_500^2) + I(NDVIBEF_2000^2)

#' Create data sets for the multivariate model
ap_cp_data <- make_multivariate_data(rhs_formula, ydata = y_train, data = train)


#THe Rhat-statistic measures the ratio of the average variance of draws within each chain to 
#the variance of the pooled draws across chains; if all chains are at equilibrium, these will 
#be the same and R^ will be one. If the chains have not converged to a common distribution, the
#R^ statistic will be greater than one (see Gelman et al. 2013, Stan Development Team 2018).
#The rule of thumb is that R-hat values for all parameters should be less than 1.1 
#Making a Rhat function
maxRhat <- function(fit) {
  max(summary(fit)$summary[, "Rhat"], na.rm = TRUE)  
}
#for rstanarm object
maxRhat_rs <- function(fit) {
  max(fit$stan_summary[, "Rhat"], na.rm = TRUE)  
}


#' univariate probit regression with rstanarm
#'Culex perexiguus
full_formula_cp <- as.formula(paste("Cxperpre ~", as.character(rhs_formula)[2]))
fit_cp <- stan_glm(full_formula_cp, data = train, refresh = 0, family = binomial(link = "probit"), 
                       init_r = 1.4, seed = 333, iter = 2000)
#Maximum Rhat-value of all parameters
maxRhat_rs(fit_cp)
#Chains seem to converge.

#Anopheles atroparvus
full_formula_at <- as.formula(paste("Anatropre ~", as.character(rhs_formula)[2]))
fit_at <- stan_glm(full_formula_at, data = train, refresh = 0, family = binomial(link = "probit"), 
                   init_r = 1.4, seed = 333, iter = 2000)
#Maximum Rhat-value of all parameters
maxRhat_rs(fit_at)
#Chains seem to converge.

#fitting the multivariate probit model
fit_mul <- stan(
  file = "multivariate_probit.stan",  # Stan program
  data = ap_cp_data,    # named list of data
  chains = 4,             # number of Markov chains
  iter = 2000,            # total number of iterations per chain
  refresh = 0,             # no progress shown
  seed = 2
)

maxRhat(fit_mul)
#Chains seem to converge

#'Prediction function for multivariate probit model with four arguments:
#'fit: fitted stanfit object
#'formula: formula object in the form ~ X1 + X2 + ...
#'newdata: Dataframe with the data specified in 'formula' (X1, X2 ... as columnnames )
#'n_max: maximum number of draws from the posterior predictive distribution

predict_mul <- function(fit, formula, newdata, n_max = 1000) {
  model_matrix <- stats::model.matrix(formula, newdata) #create appropriate matrix with all the data
  beta_dim <- fit@par_dims$beta
  omega_dim <- fit@par_dims$Omega
  stopifnot(NCOL(model_matrix) == beta_dim[2])
  stopifnot(omega_dim[1] == 2)
  #create objects with the names of the model parameters
  beta_pars <- paste0("beta[", rep(1:beta_dim[1], each = beta_dim[2]), ",",
                      rep(1:beta_dim[2], beta_dim[1]), "]")
  Omega_pars <- paste0("Omega[",  rep(1:omega_dim[1], each = omega_dim[2]), ",",
                       rep(1:omega_dim[2], omega_dim[1]), "]")
  #extract drawn parameter values from the fitted model
  draws_df <- as.data.frame(rstan::extract(fit, pars = c(beta_pars, Omega_pars)))
  #If n_max is smaller than the number of sampled parameters in the fitted model, then draw n_max-
  #times these parameters from these sampled parameters
  if (n_max < NROW(draws_df)) {
    i <- sample.int(NROW(draws_df), n_max, replace = FALSE)  
    draws_df <- draws_df[i, ]
  }
  #predict the probabilities NROW(draws-df)-times
  res <- lapply (1:NROW(draws_df), function(j) {
    #calculate the vallues of the latent variables (X %*% beta)
    means <- sapply(1:beta_dim[1], function(i) {
      betas <- as.numeric(draws_df[j, paste0("beta.", i, ".", 1:beta_dim[2], ".")])
      model_matrix %*% betas  
    })
    #obtain residual correlation parameter
    corr <- draws_df[j, "Omega.1.2."]
    #calculate the CDF for the different combinations of presence-absence (joint probabilities)
    p11 <- pbivnorm(x = means[,1], y = means[,2], rho = corr)
    p10 <- pbivnorm(x = means[,1], y = -means[,2], rho = -corr)
    p01 <- pbivnorm(x = -means[,1], y = means[,2], rho = -corr)
    p00 <- 1 - (p11 + p10 + p01)
    #calculate the different probabilities for the single species
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
##########################################In-sample (not on test data)###########################
##################################a. Comparing the coefficients: Size and credibility intervals###

#Extracting the posterior draws of the betas from the multivariate model
extract_betas <- function(fit, species = 1) {
  beta_dim <- fit@par_dims$beta
  stopifnot(species <= beta_dim[1])
  beta_pars <- paste0("beta[", species, ",", 1:beta_dim[2], "]")
  as.data.frame(rstan::extract(fit, pars = beta_pars))
}

apcp_cp <- extract_betas(fit_mul, 1)
apcp_at <- extract_betas(fit_mul, 2)

#'Plotting coefficients of uni- vs. coefficients of multivariate models and their 95%-Credibility
#'intervals: https://github.com/stan-dev/bayesplot/issues/232

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

#For Culex perexiguus
colnames(apcp_cp) <- colnames(as.matrix(fit_cp))
plot_betas(fit_cp, apcp_cp, "Univariate", "Multivariate", "Culex perexiguus")

#For Anopheles atroparvus
colnames(apcp_at) <- colnames(as.matrix(fit_at))
plot_betas(fit_at, apcp_at, "Univariate", "Multivariate", "Anopheles atroparvus")

#Mean difference of coefficients in absolute terms
#For Culex
abs(fit_cp$coefficients -  sapply(apcp_cp, median)) %>% mean
#>> Difference is very small
#For Anopheles
abs(fit_at$coefficients -  sapply(apcp_at, median)) %>% mean
#>> Difference is small

#Mean difference in relative terms (difference between every coefficient divided by average of
#the two coefficients)
#For Culex
(abs(fit_cp$coefficients -  sapply(apcp_cp, median)) / rowMeans(cbind(fit_cp$coefficients, sapply(apcp_cp, median)))) %>% abs %>% mean
#Difference is small

#For Anopheles atroparvus
(abs(fit_at$coefficients -  sapply(apcp_at, median)) / rowMeans(cbind(fit_at$coefficients, sapply(apcp_at, median)))) %>% abs %>% mean
#Difference is considerable >> Reason is the September Coefficient, But for that one the uncertainty
#also is very high >> reason: few datapoints >> So, we stick to our hypothesis that coefficients
#are the same.
#####################################6.1.b Correlation between the responses (raw vs. residual correlation)##########
#raw pearson correlation for the two species on the train set
cor(y_train$Cxperpre, y_train$Anatropre) 

#a better measure for similarity of binary data (Jaccard Index 
#(Intersection(A,B)/Union(A,B)) )
library('clusteval')
cluster_similarity(y_train$Cxperpre,y_train$Anatropre, similarity="jaccard",
                   method="independence")
#not a big difference compared to raw correlation

#Residual Correlations
#extract the posterior draws of the residual correlation (Pearson correlation coefficient)
omega <-as.data.frame(rstan::extract(fit_mul, pars = "Omega[1,2]"))
#extract the mean of the posterior draws. This is our estimate for the residual correlation
mean(omega$Omega.1.2.)
#residual correlation is smaller >> reason: shared environmental response is
#accounted for (Most coefficient have the same sign (only a couple of month
#dummies are different)! But there is still a lot of unexplained Correlation.

#####################################6.1.c. Response Curves######################################
####For Culex perexiguus
###'Response Curve für IA_500:
###'predict the probability of presence for different ia-values (over the range of ia-values in
###' the dataset), holding the other covariates at their mean (0) or at a constant factor level
###'  (Month)

#a sequence of IA_500-values from min to max with 50 steps
ia <- seq(min(df$IA_500), max(df$IA_500), length.out = 50)
#names of covariates
nc <- c("IA_500", "NDVIBEF_2000", "Month")
#the mean of the covariates is 0, bc we normalized them
#creating  x-data dataframe with all the averages (which are zero)
data <- as.data.frame(matrix(0, ncol = length(nc), nrow = length(ia)))
names(data) <- nc  
#replace IA_500 0s with the sequence from above
data$IA_500 <-ia
#Next, we want to create a dataframe with "data" for all the different months (bc we want to have
#a response curve for every month separately)
#initialize empty data frame that will be filled in the loop over the months
xdata <- as.data.frame(matrix(0, ncol = length(nc), nrow = 0))
names(xdata) <- nc  
#run a loop over all the months and append month after month (with the same 
#covariates) to the dataframe xdata
for (i in levels(df$Month)){
  d<- data
  d$Month <- i
  xdata <- rbind(xdata, d)
}
#Converting Month to a factor variable
xdata$Month <- factor(xdata$Month, levels =c("April", "May", "June", "July", "August",
                                         "September"))

#Run the posterior simulations of the univariate model with xdata, we want the predictions on
#the "probability scale". That's why we use posterior_linpred with transform = T
d_res <- posterior_linpred(fit_cp, newdata = xdata, transform = T, seed = 333)
#taking the mean of the posterior simulations as our prediction
univariate <- colMeans(d_res)


#Running the posterior simulations of the multivariate model
#Predictions for multivariate model
pred_mul_ia <- predict_mul(fit_mul, rhs_formula, xdata, n_max = 4000)

#extract the unconditional predictions for Culex and take the mean of the simulations 
#per observations as predictions
pred_ia_cp_un <- rowMeans(pred_mul_ia$p1_uncond)
#conditional predictions, if Anopheles atroparvus is present
pred_ia_cp_con1 <- rowMeans(pred_mul_ia$p1_cond)
#conditional predictions, if Anopheles atroparvus is absent
pred_ia_cp_con0 <- rowMeans(pred_mul_ia$p1_cond_0)

#####Plotting univariate and unconditional multivariate response curves in one plot
#creating a dataframe for ggplot with xdata, the univariate predictions and the unconditional
#multivariate predictions
ggd_ia <- cbind(xdata, univariate, pred_ia_cp_un)
names(ggd_ia)[length(names(ggd_ia))] <- "multivariate"

#Do the ggplot 
ggplot(data = ggd_ia, aes(x = IA_500, color = Month)) +
  geom_point(aes(y = univariate, shape = "univariate")) + 
  geom_point(aes(y = multivariate, shape = "multivariate")) + ylab("Predicted Probability") + 
  xlab("IA_500 in Standard Units") +
  ggtitle ("Univariate vs. Unconditional Multivariate Response Curves \n of Inundation Area for Culex perexiguus") +
  guides(color=guide_legend(title="Month"), shape = guide_legend(title="Model"))

#Response curves are more or less identical >> supports out hypothesis

#####Plotting conditional multivariate response curves
#combining covariates and conditional predictions
dr_con_ia <- cbind(xdata, pred_ia_cp_con1, pred_ia_cp_con0)
names(dr_con_ia)[4:5] <- c("present", "absent")

#The plot:"Response Curve of Inundation Area for Culex perexiguus in Multivariate Model Conditional on Presence of Anophles troparvus"
ggplot(data = dr_con_ia, aes(x = IA_500, color = Month)) +
  geom_point(aes(y = present, shape = "present")) + 
  geom_point(aes(y = absent, shape = "absent")) + 
  ylab("Predicted Probability") + xlab("IA_500 in Standard Units") +
  ggtitle ("Response Curve of Inundation Area for Culex perexiguus in Multivariate Model \n Conditional on Presence of Anophles troparvus") +
  guides(color=guide_legend(title="Month"), shape = guide_legend(title = "Anopheles troparvus"))

#interesting results, ordering of the effect of month on pred, changes depending on whether AT is present
#'or not. Does that make sense? Yes, because P(A|B) changes if mean(latent variable z)/
#'mean(P(B)) changes and this can change with the months. 
#'
####' Plot that shows that average effect of presence of anopheles per month depends
####' on the mean presence of Anopheles in that month
#Make a dataframe with all the continuous predictors at 0 for each month
d_dif <- as.data.frame(matrix(0, ncol = length(nc), nrow = length(levels(df$Month))))
names(d_dif) <- nc
d_dif$Month <- levels(df$Month)
#predict probabilities for multivariate model
preds <- predict_mul(fit_mul, rhs_formula, d_dif, n_max = 4000)

#calculate difference between predictions with present anopheles and absent anopheles per month
dif <- rowMeans(preds$p1_cond) - rowMeans(preds$p1_cond_0)

#make dataframe for plotting with this difference, the month and the average presence per month
#in the train set
d_dif <- data.frame(levels(df$Month), dif, tapply(df$Anatropre, df$Month, mean))
names(d_dif) <- c("Month", "difference", "MeanPresence")
#make the plot
ggplot(data = d_dif, aes(x = MeanPresence, y = difference, color = Month)) + 
  geom_point() + labs(title = "Mean Presence of Anopheles atroparvus per Month vs. Difference in Predicted Probability of \nCulex perexiguus between Observations with Present and Absent Anopheles atroparvus")
#At least, for september our suspicion holds true that a low mean presence is associated with 
#a high condtional probability of Culex given the presence of Anopheles. For the other months
#there is no trend detectable between difference and MeanPresence.

####Doing a plot of multivariate conditional predictions and univariate prediction per month

#make a list of plots with one plot for every month. Each plot contains response three response
#curves: (i) univariate model, (ii) conditional predictions with absent Anopheles and (iii)
#conditional predictions with present Anopheles
#Preparing the input data
plots_uv_con <- vector(mode = "list", length = length(levels(xdata$Month)))
for( i in levels(xdata$Month)){
  #index variable
  j <-match(i, levels(xdata$Month))
  #Number of observations
  nu <- nrow(xdata[xdata["Month"] == i,])
  #three-times the xdata in one data frame, so that every prediction type is represented
  d <- rbind(xdata[xdata["Month"] == i,], xdata[xdata["Month"] == i,], xdata[xdata["Month"] == i,])
  #creating one column that has the presence-absence of Anopheles for the multivariate models and for the univariate
  #model the variable takes on another factor level (2)
  d$mode <- factor(c(rep(0, nu), rep(1, nu), rep(2, nu)), labels = c("Multivariate Model with absent Anopheles", 
                                                                     "Multivariate Model with present Anopheles",
                                                                     "Univariate Model"))
  #adding the predictions
  d$pred <- c(dr_con_ia$absent[dr_con_ia$Month == i], dr_con_ia$present[dr_con_ia$Month == i],
              ggd_ia$univariate[ggd_ia$Month == i])
  
  #Doing the ggplot
  plot <- ggplot(data = d, aes(x = IA_500, y = pred, shape = mode)) +
    geom_point() + ylab("Predicted Probability") + xlab("IA_500 in Standard Units") +
    ggtitle (paste0("Response Curve of Inundation Area for Culex perexiguus in Multivariate Model \n Conditional on Presence of Anophles troparvus vs. Univariate Model for Month ", i)) +
    theme(legend.title = element_blank())
  plots_uv_con[[j]] <- plot
}
#Plot for June
plots_uv_con[[3]]
#make sense

########for NDVIBEF

####for univariate predictions
#a sequence of NDVIBEF-values from min to max with 50 steps
nd <- seq(min(df$NDVIBEF_2000), max(df$NDVIBEF_2000), length.out = 50)

#the mean of the covariates is 0, bc I normalized them
#creating  x-data dataframe with all the averages (which are zero)
data <- as.data.frame(matrix(0, ncol = length(nc), nrow = length(nd)))
names(data) <- nc  
#replace NDVIBEF_2000 0s with the sequence from above
data$NDVIBEF_2000 <-nd
#Next, we want to create a dataframe with "data" for all the different months (bc we want to have
#a response curve for every month separately)
#initialize empty data frame that will be filled in the loop over the months
xdata <- as.data.frame(matrix(, ncol = length(nc), nrow = 0))
names(xdata) <- nc  
#run a loop over all the months and append month after month (with the same 
#covariates) to the dataframe xdata
for (i in levels(df$Month)){
  d<- data
  d$Month <- i
  xdata <- rbind(xdata, d)
}
#Converting Month to a factor variable
xdata$Month <- factor(xdata$Month, levels =c("April", "May", "June", "July", "August",
                                         "September"))
#Run the posterior simulations with xdata
d_res <- posterior_linpred(fit_cp, newdata = xdata, seed = 333, transform = T)
#getting the mean predictions as our predicted probability
univariate <- colMeans(d_res)

#Running the posterior simulations of the multivariate model
#Predictions for multivariate model
pred_mul_nd <- predict_mul(fit_mul, rhs_formula, xdata, n_max = 4000)

#extract the unconditional predictions for Culex and take the mean of the simulations 
#per observations as predictions
pred_nd_cp_un <- rowMeans(pred_mul_nd$p1_uncond)
#conditional predictions, if Anopheles atroparvus is present
pred_nd_cp_con1 <- rowMeans(pred_mul_nd$p1_cond)
#conditional predictions, if Anopheles atroparvus is absent
pred_nd_cp_con0 <- rowMeans(pred_mul_nd$p1_cond_0)

#####Plotting univariate and unconditional multivariate response curves in one plot
#creating a dataframe for ggplot with xdata, the univariate predictions and the unconditional
#multivariate predictions
ggd_nd <- cbind(xdata, univariate, pred_nd_cp_un)
names(ggd_nd)[length(names(ggd_nd))] <- "multivariate"

#Do the ggplot 
ggplot(data = ggd_nd, aes(x = NDVIBEF_2000, color = Month)) +
  geom_point(aes(y = univariate, shape = "univariate")) + 
  geom_point(aes(y = multivariate, shape = "multivariate")) + ylab("Predicted Probability") + 
  xlab("NDVIBEF_2000 in Standard Units") +
  ggtitle ("Univariate vs. Unconditional Multivariate Response Curves \n of \"NDVI one Month Before\"  for Culex perexiguus") +
  guides(color=guide_legend(title="Month"), shape = guide_legend(title="Model"))

#'Response curves are more or less identical >> supports out hypothesis! However, the overall 
#'shape is a little weird. Contradicts the niche hypothesis right? But remember it's NDVI >> No
#'linear interpretation I guess.

#####Plotting conditional multivariate response curves
#combining covariates and conditional predictions
dr_con_nd <- cbind(xdata, pred_nd_cp_con1, pred_nd_cp_con0)
names(dr_con_nd)[4:5] <- c("present", "absent")

#Plotting
ggplot(data = dr_con_nd, aes(x = NDVIBEF_2000, color = Month)) +
  geom_point(aes(y = present, shape = "present")) + 
  geom_point(aes(y = absent, shape = "absent")) + 
  ylab("Predicted Probability") + xlab("NDVIBEF_2000 in Standard Units") +
  ggtitle ("Multivariate Response Curve of \" NDVI Month before \" for Culex perexiguus  \n Conditional on Presence of Anophles troparvus") +
  guides(color=guide_legend(title="Month"), shape = guide_legend(title = "Anopheles troparvus"))

#interesting results, ordering of the effect of month on pred, changes depending on whether AT is present
#'or not. 

####Doing a plot of multivariate conditional predictions and univariate prediction per month

#make a list of plots with one plot for every month. Each plot contains response three response
#curves: (i) univariate model, (ii) conditional predictions with absent Anopheles and (iii)
#conditional predictions with present Anopheles
#Preparing the input data
plots_uv_con <- vector(mode = "list", length = length(levels(xdata$Month)))
for( i in levels(xdata$Month)){
  #index variable
  j <-match(i, levels(xdata$Month))
  #Number of observations
  nu <- nrow(xdata[xdata["Month"] == i,])
  #three-times the xdata in one data frame, so that every prediction type is represented
  d <- rbind(xdata[xdata["Month"] == i,], xdata[xdata["Month"] == i,], xdata[xdata["Month"] == i,])
  #creating one column that has the presence-absence of Anopheles for the multivaraite models and for the univariate
  #model the variable takes on another factor level (2)
  d$mode <- factor(c(rep(0, nu), rep(1, nu), rep(2, nu)), labels = c("Multivariate Model with absent Anopheles", 
                                                                     "Multivariate Model with present Anopheles",
                                                                     "Univariate Model"))
  #adding the predictions
  d$pred <- c(dr_con_nd$absent[dr_con_nd$Month == i], dr_con_nd$present[dr_con_nd$Month == i],
              ggd_nd$univariate[ggd_nd$Month == i])
  
  #Doing the ggplot
  plot <- ggplot(data = d, aes(x = NDVIBEF_2000, y = pred, shape = mode)) +
    geom_point() + ylab("Predicted Probability") + xlab("NDVIBEF_2000 in Standard Units") +
    ggtitle (paste0("Response Curve of \" NDVI Month Before \" for Culex perexiguus in Multivariate Model \n Conditional on Presence of Anophles troparvus vs. Univariate Model for Month ", i)) +
    theme(legend.title = element_blank())
  plots_uv_con[[j]] <- plot
}
#Plot for June
plots_uv_con[[5]]
#'For extreme positive values, univariate predictions approach multivariate predictions with
#'present Anopheles

####For Anopheles atroparvus ####
###'Response Curve für IA_500:
###'predict the probability of presence for different ia-values (over the range of ia-values in
###' the dataset), holding the other covariates at their mean (0) or at a constant factor level
###'  (Month)

#creating  x-data dataframe with all the averages (which are zero)
data <- as.data.frame(matrix(0, ncol = length(nc), nrow = length(ia)))
names(data) <- nc  
#replace IA_500 0s with the sequence from above
data$IA_500 <-ia
#Next, we want to create a dataframe with "data" for all the different months (bc we want to have
#a response curve for every month separately)
#initialize empty data frame that will be filled in the loop over the months
xdata <- as.data.frame(matrix(0, ncol = length(nc), nrow = 0))
names(xdata) <- nc  
#run a loop over all the months and append month after month (with the same 
#covariates) to the dataframe xdata
for (i in levels(df$Month)){
  d<- data
  d$Month <- i
  xdata <- rbind(xdata, d)
}
#Converting Month to a factor variable
xdata$Month <- factor(xdata$Month, levels =c("April", "May", "June", "July", "August",
                                             "September"))

#Run the posterior simulations of the univariate model with xdata, we want the predictions on
#the "probability scale". That's why we use posterior_linpred with transform = T
d_res <- posterior_linpred(fit_at, newdata = xdata, transform = T, seed = 333)
#taking the mean of the posterior simulations as our prediction
univariate <- colMeans(d_res)


#extract the unconditional predictions for Culex and take the mean of the simulations 
#per observations as predictions
pred_ia_at_un <- rowMeans(pred_mul_ia$p2_uncond)
#conditional predictions, if Culex perexiguus is present
pred_ia_at_con1 <- rowMeans(pred_mul_ia$p2_cond)
#conditional predictions, if Culex perexiguus is absent
pred_ia_at_con0 <- rowMeans(pred_mul_ia$p2_cond_0)

#####Plotting univariate and unconditional multivariate response curves in one plot
#creating a dataframe for ggplot with xdata, the univariate predictions and the unconditional
#multivariate predictions
ggd_ia <- cbind(xdata, univariate, pred_ia_at_un)
names(ggd_ia)[length(names(ggd_ia))] <- "multivariate"

#Do the ggplot 
ggplot(data = ggd_ia, aes(x = IA_500, color = Month)) +
  geom_point(aes(y = univariate, shape = "univariate")) + 
  geom_point(aes(y = multivariate, shape = "multivariate")) + ylab("Predicted Probability") + 
  xlab("IA_500 in Standard Units") +
  ggtitle ("Univariate vs. Unconditional Multivariate Response Curves \n of Inundation Area for Anopheles atroparvus") +
  guides(color=guide_legend(title="Month"), shape = guide_legend(title="Model"))

#Response curves are more or less identical >> supports out hypothesis

#####Plotting conditional multivariate response curves
#combining covariates and conditional predictions
dr_con_ia <- cbind(xdata, pred_ia_at_con1, pred_ia_at_con0)
names(dr_con_ia)[4:5] <- c("present", "absent")

#The plot:"Response Curve of Inundation Area for Anopheles atroparvus in Multivariate Model Conditional on Presence of Culex"
ggplot(data = dr_con_ia, aes(x = IA_500, color = Month)) +
  geom_point(aes(y = present, shape = "present")) + 
  geom_point(aes(y = absent, shape = "absent")) + 
  ylab("Predicted Probability") + xlab("IA_500 in Standard Units") +
  ggtitle ("Response Curve of Inundation Area for Anopheles atroparvus in Multivariate Model \n Conditional on Presence of Culex perexiguus") +
  guides(color=guide_legend(title="Month"), shape = guide_legend(title = "Anopheles troparvus"))

#results are as expected 

####Doing a plot of multivariate conditional predictions and univariate prediction per month

#make a list of plots with one plot for every month. Each plot contains response three response
#curves: (i) univariate model, (ii) conditional predictions with absent Anopheles and (iii)
#conditional predictions with present Anopheles
#Preparing the input data
plots_uv_con <- vector(mode = "list", length = length(levels(xdata$Month)))
for( i in levels(xdata$Month)){
  #index variable
  j <-match(i, levels(xdata$Month))
  #Number of observations
  nu <- nrow(xdata[xdata["Month"] == i,])
  #three-times the xdata in one data frame, so that every prediction type is represented
  d <- rbind(xdata[xdata["Month"] == i,], xdata[xdata["Month"] == i,], xdata[xdata["Month"] == i,])
  #creating one column that has the presence-absence of Anopheles for the multivariate models and for the univariate
  #model the variable takes on another factor level (2)
  d$mode <- factor(c(rep(0, nu), rep(1, nu), rep(2, nu)), labels = c("Multivariate Model with absent Anopheles", 
                                                                     "Multivariate Model with present Anopheles",
                                                                     "Univariate Model"))
  #adding the predictions
  d$pred <- c(dr_con_ia$absent[dr_con_ia$Month == i], dr_con_ia$present[dr_con_ia$Month == i],
              ggd_ia$univariate[ggd_ia$Month == i])
  
  #Doing the ggplot
  plot <- ggplot(data = d, aes(x = IA_500, y = pred, shape = mode)) +
    geom_point() + ylab("Predicted Probability") + xlab("IA_500 in Standard Units") +
    ggtitle (paste0("Response Curve of Inundation Area for Anopheles atroparvus in Multivariate Model \n Conditional on Presence of Anophles troparvus vs. Univariate Model for Month ", i)) +
    theme(legend.title = element_blank())
  plots_uv_con[[j]] <- plot
}
#Plot for June
plots_uv_con[[3]]
#make sense

########for NDVIBEF

####for univariate predictions
#the mean of the covariates is 0, bc I normalized them
#creating  x-data dataframe with all the averages (which are zero)
data <- as.data.frame(matrix(0, ncol = length(nc), nrow = length(nd)))
names(data) <- nc  
#replace NDVIBEF_2000 0s with the sequence from above
data$NDVIBEF_2000 <-nd
#Next, we want to create a dataframe with "data" for all the different months (bc we want to have
#a response curve for every month separately)
#initialize empty data frame that will be filled in the loop over the months
xdata <- as.data.frame(matrix(, ncol = length(nc), nrow = 0))
names(xdata) <- nc  
#run a loop over all the months and append month after month (with the same 
#covariates) to the dataframe xdata
for (i in levels(df$Month)){
  d<- data
  d$Month <- i
  xdata <- rbind(xdata, d)
}
#Converting Month to a factor variable
xdata$Month <- factor(xdata$Month, levels =c("April", "May", "June", "July", "August",
                                             "September"))
#Run the posterior simulations with xdata
d_res <- posterior_linpred(fit_at, newdata = xdata, seed = 333, transform = T)
#getting the mean predictions as our predicted probability
univariate <- colMeans(d_res)

#extract the unconditional predictions for Anopheles and take the mean of the simulations 
#per observations as predictions
pred_nd_at_un <- rowMeans(pred_mul_nd$p2_uncond)
#conditional predictions, if Culex perexiguus is present
pred_nd_at_con1 <- rowMeans(pred_mul_nd$p2_cond)
#conditional predictions, if Culex perexiguus is absent
pred_nd_at_con0 <- rowMeans(pred_mul_nd$p2_cond_0)

#####Plotting univariate and unconditional multivariate response curves in one plot
#creating a dataframe for ggplot with xdata, the univariate predictions and the unconditional
#multivariate predictions
ggd_nd <- cbind(xdata, univariate, pred_nd_at_un)
names(ggd_nd)[length(names(ggd_nd))] <- "multivariate"

#Do the ggplot 
ggplot(data = ggd_nd, aes(x = NDVIBEF_2000, color = Month)) +
  geom_point(aes(y = univariate, shape = "univariate")) + 
  geom_point(aes(y = multivariate, shape = "multivariate")) + ylab("Predicted Probability") + 
  xlab("NDVIBEF_2000 in Standard Units") +
  ggtitle ("Univariate vs. Unconditional Multivariate Response Curves \n of \"NDVI one Month Before\"  for Anopheles atroparvus") +
  guides(color=guide_legend(title="Month"), shape = guide_legend(title="Model"))

#'Response curves are more or less identical >> supports out hypothesis!

#####Plotting conditional multivariate response curves
#combining covariates and conditional predictions
dr_con_nd <- cbind(xdata, pred_nd_at_con1, pred_nd_at_con0)
names(dr_con_nd)[4:5] <- c("present", "absent")

#Plotting
ggplot(data = dr_con_nd, aes(x = NDVIBEF_2000, color = Month)) +
  geom_point(aes(y = present, shape = "present")) + 
  geom_point(aes(y = absent, shape = "absent")) + 
  ylab("Predicted Probability") + xlab("NDVIBEF_2000 in Standard Units") +
  ggtitle ("Multivariate Response Curve of \" NDVI Month before \" for Anopheles atroparvus  \n Conditional on Presence of Culex perexiguus") +
  guides(color=guide_legend(title="Month"), shape = guide_legend(title = "Culex perexiguus"))

#results as expected

####Doing a plot of multivariate conditional predictions and univariate prediction per month

#make a list of plots with one plot for every month. Each plot contains response three response
#curves: (i) univariate model, (ii) conditional predictions with absent Anopheles and (iii)
#conditional predictions with present Anopheles
#Preparing the input data
plots_uv_con <- vector(mode = "list", length = length(levels(xdata$Month)))
for( i in levels(xdata$Month)){
  #index variable
  j <-match(i, levels(xdata$Month))
  #Number of observations
  nu <- nrow(xdata[xdata["Month"] == i,])
  #three-times the xdata in one data frame, so that every prediction type is represented
  d <- rbind(xdata[xdata["Month"] == i,], xdata[xdata["Month"] == i,], xdata[xdata["Month"] == i,])
  #creating one column that has the presence-absence of Culex for the multivaraite models and for the univariate
  #model the variable takes on another factor level (2)
  d$mode <- factor(c(rep(0, nu), rep(1, nu), rep(2, nu)), labels = c("Multivariate Model with absent Culex", 
                                                                     "Multivariate Model with present Culex",
                                                                     "Univariate Model"))
  #adding the predictions
  d$pred <- c(dr_con_nd$absent[dr_con_nd$Month == i], dr_con_nd$present[dr_con_nd$Month == i],
              ggd_nd$univariate[ggd_nd$Month == i])
  
  #Doing the ggplot
  plot <- ggplot(data = d, aes(x = NDVIBEF_2000, y = pred, shape = mode)) +
    geom_point() + ylab("Predicted Probability") + xlab("NDVIBEF_2000 in Standard Units") +
    ggtitle (paste0("Response Curve of \" NDVI Month Before \" for Anopheles atroparvus in Multivariate Model \n Conditional on Presence of Culex perexiguus vs. Univariate Model for Month ", i)) +
    theme(legend.title = element_blank())
  plots_uv_con[[j]] <- plot
}
#Plot for June
plots_uv_con[[3]]
#'For extreme positive values, univariate predictions approach multivariate predictions with
#'present Culex

#####################################6.1.d. Variable Importance######
####For Culex perexiguus
####for univariate model

#Define our covariates as xdata
xdata <- train[,c("Month", "IA_500", "NDVIBEF_2000")]

# create custom predict function for rstanarm "posterior_linpred" 
#We take the mean of the simulations as the prediction on the probability-scale.
pred_uv <- function(model, newdata)  {
  return(posterior_linpred(model, newdata, transform = T, seed = 333) %>% apply(2, mean))
}
#create the explain object (core object of DALEX) which contains the model, the data and
#the predict function
dal_cp <- explain(fit_cp, xdata, y = train$Cxperpre, predict_function = pred_uv,
                  type = "classification", label = "Univariate Probit")

# calculate the permutation-based variable importance measure (Here: the difference in
#(1-AUC) between original data and permuted data per covariate)
set.seed(1980)
vi_cp <- model_parts(dal_cp, type = "difference", B = 50)

#plot the results
plot(vi_cp) + 
  labs(title = "Variable Importance over 50 Permuations", subtitle = "created for the univariate probit model of Culex perexiguus") 
#So, Month is the  most important variable, followed by IA_500 and at last 
#NDVIBEF_2000. When permuting the Mes variable entries and then predicting
#our response, the resulting AUC is roughly .11 worse than for the 
#predictions without the permutations.

####For multivariate unconditional predictions

# create custom predict function:
pred_mv_cp <- function(model, newdata)  {
  #run the predict function
  pred <- predict_mul(model, rhs_formula, newdata, n_max = 4000)
  #return the means of the simulations as the unconditional predictions of species 1 (Culex perexiguus)
  return(rowMeans(pred$p1_uncond))
}

#make the explainer
dal_cpmv <- explain(fit_mul, xdata, y = train$Cxperpre,
                    predict_function = pred_mv_cp, type = "classification",
                    label = "Unconditional Multivariate Probit")

#permutation-based variable importance measure
set.seed(1980)
vi_cpmv <- model_parts(dal_cpmv, type = "difference", B = 50) 

#plot the results
plot(vi_cpmv) + labs(title = "Variable Importance over 50 Permutations", subtitle = "created for the multivariate probit model of Culex perexiguus")
#Most important variable is Month, then IA_500 and then NDVIBEF_2000. Very similiar to univariate
#case

####plotting the two: 
#making a tibble with (i) a variable column indicating which variable is permuted, (ii) a loss column giving the 
#(1-Auc)-losses and (iii) a type column indicating the model type
d_gg <- tibble(var = c(vi_cp$variable, vi_cpmv$variable), loss = c(vi_cp$dropout_loss, vi_cpmv$dropout_loss), type = c(vi_cp$label, vi_cpmv$label))
#converting var and type to a factor
d_gg$var <- factor(d_gg$var, levels = c("NDVIBEF_2000", "IA_500", "Month", "_baseline_", "_full_model_"))
d_gg$type <- factor(d_gg$type, levels = c("Univariate Probit", "Unconditional Multivariate Probit"))

#deleting baseline and full model (We do not need these loss values)
d_gg <-d_gg[!(d_gg$var == "_full_model_" | d_gg$var == "_baseline_"),]

#set overall theme of ggplotting to grey
theme_set(theme_grey())
#Plot "Variable Importance Boxplots"
ggplot(d_gg, aes(x = loss, y = var)) +
  geom_boxplot(aes(color = type)) + # Boxplot shows Boxes: 25% and 75 % Quantile; vertical line: median; lower whisker = smallest observation greater than or equal to lower hinge - 1.5 * IQR
  facet_wrap(~type, nrow =3) + xlim(0,.2) +
  labs(title = "Variable Importance Boxplots",
       subtitle = "Created for multi- and univariate models of Culex perexiguus with 50 permutation runs", y = "Variable",
       x = "(1 - AUC)-Loss after Permutations") +
  # Suppress the legend since color isn't actually providing any information
  theme(legend.position = "none") 

#### Conditional multivariate predictions

#xdata with PA of Anopheles. This way, we treat the PA-data like a variable for which we can
#calculate VI
xdata_con <- as_tibble(cbind(xdata, y_train[,2]))

#make the predict function, including PA-anopheles as a covariate to
#permute/manipulate
pred_con <- function(model, newdata){
  #run the predict function
  pred <- predict_mul(model, rhs_formula, newdata, n_max = 4000)
  #'return the means of the simulations for the conditional predictions with present Anopheles
  #'if Anopheles is present, return the mean of the simulations with absent Anopheles, if
  #'Anopheles is absent
  return(ifelse(newdata$Anatropre == 1, rowMeans(pred$p1_cond), rowMeans(pred$p1_cond_0)))
}

#make the explainer
dal_cp_con <- explain(fit_mul, xdata_con, y = train$Cxperpre,
                      predict_function = pred_con, type = "classification")

#permutation-based variable importance measure
set.seed(1980)
#Doing it for 50 permutations
vi_cp_con <- model_parts(dal_cp_con, type = "difference", B = 50)

#plotting the results
plot(vi_cp_con) + labs(title = "Variable Importance", subtitle = "created for the multivariate probit model of Culex perexiguus conditional on Anopheles troparvus")

#'Anopheles is the most important variable, followed by IA_500 and Month. So, month and IA_500
#'swap places compared to the unconditional case. Their magnitude is way smaller.
#The reason is that the PA-data takes away some of the predictive power of the environmental 
#covariates and especially of the month variable.

####For Anopheles####
###For Univariate Probit Model

#make the explainer
dal_at <- explain(fit_at, xdata, y = train$Anatropre, predict_function = pred_uv,
                  type = "classification", label = "Univariate Probit")


# calculate the permutation-based variable importance measure 
set.seed(1980)
vi_at <- model_parts(dal_at, type = "difference", B = 50)
#plot it
plot(vi_at) + 
  labs(title = "Variable Importance over 50 Permuations", subtitle = "created for the univariate probit model of Anopheles troparvus") 
#So, Month is the  most important variable, followed by NDVIBEF_2000 and at last 
#IA_500. When permuting the Mes variable entries and then predicting
#our response, the resulting AUC is roughly .2 worse than for the 
#predictions without the permutations. 

####For unconditional predictions in multivariate model

#We need to slightly change our predict function (just the last line) compared to the Culex case
pred_mv_at <- function(model, newdata)  {
  #run the predict function
  pred <- predict_mul(model, rhs_formula, newdata, n_max = 4000)
  #return the means of the simulation for the unconditional predictions of species 2 (Anopheles atroparvus)
  return(rowMeans(pred$p2_uncond))
}

#make the explainer
dal_mvat <- explain(fit_mul, xdata, y = train$Anatropre,
                    predict_function = pred_mv_at, type = "classification", 
                    label = "Unconditional Multivariate Probit")

#permutation-based variable importance measure
set.seed(1980)
vi_mvat <- model_parts(dal_mvat, type = "difference", B = 50) 
#pretty similar to univariate model!

#plot the results
plot(vi_mvat) + labs(title = "Variable Importance over 50 Permutations", subtitle = "created for the multivariate probit model of Anopheles troparvus")
#Most important variable is Mes, then NDVI_BEF and then IA, similar results as for the univariate
#model

###Plotting VI of both models in one plot

# making a dataframe with all the important data: 1. The variable which is permuted, 2. the drop 
#out loss for every run and 3. the according model type labels
d_ggat <- tibble(var = c(vi_at$variable, vi_mvat$variable),
                 loss = c(vi_at$dropout_loss, vi_mvat$dropout_loss),
                 type = c(vi_at$label, vi_mvat$label))
#converting var and type to a factor
d_ggat$var <- factor(d_ggat$var, levels = c("IA_500", "NDVIBEF_2000", "Month", "_full_model_", 
                                            "_baseline_"))
d_ggat$type <- factor(d_ggat$type, levels = c("Univariate Probit", "Unconditional Multivariate Probit"))

#deleting baseline and full model 
d_ggat <-d_ggat[!(d_ggat$var == "_full_model_" | d_ggat$var == "_baseline_"),]

#Plot the comparison
ggplot(d_ggat, aes(x = loss, y = var)) +
  geom_boxplot(aes(color = type)) + # Boxplot shows Boxes: 25% and 75 % Quantile; vertical line: median; lower whisker = smallest observation greater than or equal to lower hinge - 1.5 * IQR
  facet_wrap(~type, nrow =3) +
  labs(title = "Variable Importance Boxplots",
       subtitle = "Created for multi- and univariate models of Anopheles troparvus with 50 permutation runs", y = "Variable",
       x = "(1 - AUC)-Loss after Permutations") +
  # Suppress the legend since color isn't actually providing any information
  theme(legend.position = "none") 
#look very similiar
##################Doing it for conditional predictions
#Treat the conditioning as a covariate (inlcude it as one variable in your
#variable importance analysis)

#xdata with PA of Culex
xdata_con <- as_tibble(cbind(xdata, y_train[,1]))

#make the predict function, including PA-Culex as a covariate to
#permute/manipulate
pred_atcon <- function(model, newdata){
  #run the predict function
  pred <- predict_mul(model, rhs_formula, newdata, n_max = 4000)
  #'return the means of the simulations for the conditional predictions with present Culex
  #'if Culex is present, return the mean of the simulations with absent Culex, if
  #'Culex is absent
  return(ifelse(newdata$Cxperpre == 1, rowMeans(pred$p2_cond), rowMeans(pred$p2_cond_0)))
}

#make the explainer
dal_at_con <- explain(fit_mul, xdata_con, y = train$Anatropre,
                      predict_function = pred_atcon, type = "classification")

#permutation-based variable importance measure
set.seed(1980)
#Doing it for 50 permutations
vi_at_con <- model_parts(dal_at_con, type = "difference", B = 50)

#plotting the results
plot(vi_at_con) + labs(title = "Variable Importance", subtitle = "created for the multivariate probit model of Anopheles troparvus conditional on Culex perexiguus")
#Culex isnt that important; is that in accordance with the AUCs? AUC is .04
#better for conditional model. So, here VI of conditioning lower way lower as in the case above,
#but later the improvement of AUC on test set in absolute terms is the same.
#Mes variable does not lose much of its importance measured in dif. The other two variables lose 
#considerably.
##########################################6.2. Out-of-Sample##########################
#####################################6.2.a. Conditional Predictions of multivariate model vs. Unconditional Predictions of multivariate model vs. predictions of univariate model####

#We evaluate the predictive performance on our test set with the AUC.

#####For Culex perexiguus
#We make the three predictions:(i) univariate, (ii) unconditional multivariate and (iii)
#conditional multivariate

###(i) univariate predictions on the probability scale
pred_cp_uv <- posterior_linpred(fit_cp, newdata = test, transform = T, seed = 333)
#take the average of these draws per observation as the predicted probability
pred_cp_uv <- colMeans((pred_cp_uv))
#AUC
auc_cp_uv <- auc(response = test$Cxperpre, predictor = pred_cp_uv)
auc_cp_uv
#a AUC of .76 >> is ok for our purposes

###(ii) Unconditional multivariate
#predictions
pred_mv_test <- predict_mul(fit_mul, rhs_formula, test, n_max = 4000)
pred_cp_un <- rowMeans(pred_mv_test$p1_uncond)
#AUC
auc_cp_mvun <- auc(response = test$Cxperpre, predictor = pred_cp_un)
auc_cp_mvun
#AUC roughly .76 >> the same as the univariate model >> corroborates our hypothesis

###Plot the univariate conditions against the unconditional multivariate ones
d_gg_uc_cp <- data.frame(cbind(pred_cp_uv, pred_cp_un, test$Cxperpre))
names(d_gg_uc_cp) <- c("cp_uv", "cp_mv_un", "cp")
#make a factor out of cp. This comes in handy for the plotting
d_gg_uc_cp$cp <- factor(d_gg_uc_cp$cp, levels = c(0, 1), labels = c("absent", "present"))

#Plot "Univariate vs. Unconditional Multivariate Predictions for Culex perexiguus"
ggplot(d_gg_uc_cp, aes(x = cp_uv, y = cp_mv_un, color = cp)) + geom_point() +
  ggtitle("Univariate vs. Unconditional Multivariate Predictions for Culex perexiguus") +
  xlab("Predictions from Univariate Probit ") + geom_abline(slope=1, intercept=0) + 
  ylab("Unconditional Predictions from Multivariate Probit") +
  labs(color = "Observed PA of Culex perexiguus") + ylim(0, 1) 

#Points are centered around the identity line, as we would expect according to hypothesis 3.

###(iii)Conditional predictions
#Culex perexiguus conditioned on Anopheles troparvus

#retrieving the conditional predictions
p_cp_mvco <- ifelse(test$Anatropre == 1, rowMeans(pred_mv_test$p1_cond), rowMeans(pred_mv_test$p1_cond_0))
#AUC
auc_cp_mvco <- auc(response = test$Cxperpre, predictor = p_cp_mvco)
auc_cp_mvco
#AUC roughly .8 >> So, there is a considerable improvement of .04 in AUC >> corroborates
#our hypothesis that PA data enhances predictions :)

###Plotting Univariate Predictions vs. Conditional Multivariate Predictions
d_gg_cp <- data.frame(cbind(pred_cp_uv, p_cp_mvco, test$Anatropre, test$Cxperpre))
names(d_gg_cp) <- c("cp_uv", "cp_mv", "at", "cp")
#Make the PA-variables to factors
d_gg_cp$cp <- factor(d_gg_cp$cp, levels = c(0, 1), labels = c("absent", "present"))
d_gg_cp$at <- factor(d_gg_cp$at, levels = c(0, 1), labels = c("absent", "present"))

ggplot(d_gg_cp, aes(x=cp_uv, y=cp_mv, color=at, shape = cp)) + geom_point() +
  ggtitle("Univariate vs. Conditional Multivariate Predictions for Culex perexiguus") +
  xlab("Predictions from Univariate Probit ") + ylab("Conditional Predictions from Multivariate Probit") +
  labs( color = "PA of Anopheles troparvus", shape = "PA of Culex perexiguus") + ylim(0,1) +
  geom_abline(slope=1, intercept=0)

#You can see that there is a positive linear relationship between the predictions
#of the two models grouped by the PA of Anopheles troparvus (the species we conditioned on).
#This indicates that both models roughly do the same/environmental signals are treated 
#similarily. The slopes of the two lines are a little less than  1 (we would expect 1, right?).
#We can still the very slight effect that the univariate models tends to predict extremer values.
#Or is this just by chance?
#You can see that the conditioning on Anopheles troparvus has a 
#clear effect (The blue and red points form two distinct groups) on the multivariate 
#predictions compared to the univariate predictions. The multivariate model predicts a roughly
#12.5 % points higher probability of occurence for plots where Anopheles troparvus is present 
#compared to plots where it's absent.

####for Anopheles ####
#We evaluate predictive performance on our test set with the AUC
#We make the three predictions ((i) univariate, (ii) unconditional multivariate and (iii)
#conditional multivariate)

###(i) univariate predictions on the probability scale
pred_at_uv <- posterior_linpred(fit_at, newdata = test, transform = T, seed = 333)
#take the mean of the simulation and take it as the prediction for the particular observation
pred_at_uv <- colMeans(pred_at_uv)
#AUC
auc_at_uv <- auc(response = test$Anatropre, predictor = pred_at_uv)
auc_at_uv
#a AUC of .85 >> is ok for our purposes

###Unconditional multivariate predictions

#retrieving the unconditional predictions (the mean of the simulations)
p_at_mv_un <- rowMeans(pred_mv_test$p2_uncond)
#AUC
auc_at_mvun <- auc(response = test$Anatropre, predictor = p_at_mv_un)
auc_at_mvun
#AUC roughly .85 >> roughly the same as the univariate model >> corroborates our hypothesis

###Plot the univariate against the unconditional multivariate predictions
#Preparing the data
d_gg_uc_at <- data.frame(cbind(pred_at_uv, p_at_mv_un, y_test$Anatropre))
names(d_gg_uc_at) <- c("at_uv", "at_mv_un", "at")
d_gg_uc_at$at <- factor(d_gg_uc_at$at, levels = c(0, 1), labels = c("absent", "present"))
#Plotting "Univariate vs. Unconditional Multivariate Predictions for Anopheles troparvus"
ggplot(d_gg_uc_at, aes(x=at_uv, y=at_mv_un, color=factor(at))) + geom_point() +
  ggtitle("Univariate vs. Unconditional Multivariate Predictions for Anopheles troparvus") +
  xlab("Predictions from Univariate Probit ") + 
  ylab("Unconditional Predictions from Multivariate Probit") +
  labs(color = "True PA of Culex perexiguus") + ylim(0, 1) + geom_abline(slope=1, intercept=0)

#Predictions are on the identity line >> corroborates hypothesis 3

###Conditional predictions
#Anopheles troparvus conditioned on Culex perexiguus 
#retrieving the conditional predictions
p_at_mvco <- ifelse(test$Cxperpre == 1, rowMeans(pred_mv_test$p2_cond), rowMeans(pred_mv_test$p2_cond_0))

#AUC
auc_at_mvco <- auc(response = test$Anatropre, predictor = p_at_mvco)
auc_at_mvco
#AUC roughly .89 >> So, there is a considerable improvement of .04 in AUC >> corroborates
#our hypothesis that PA-data of other species enhances predictions.

#Plotting Univariate Predictions vs. Conditional Multivariate Predictions
d_gg_at <- data.frame(cbind(pred_at_uv, p_at_mvco, test$Anatropre, test$Cxperpre))
names(d_gg_at) <- c("at_uv", "at_mv", "at", "cp")
#Make the PA-variables to factors
d_gg_at$cp <- factor(d_gg_at$cp, levels = c(0, 1), labels = c("absent", "present"))
d_gg_at$at <- factor(d_gg_at$at, levels = c(0, 1), labels = c("absent", "present"))

ggplot(d_gg_at, aes(x=at_uv, y=at_mv, color=cp, shape = at)) + geom_point() +
  ggtitle("Univariate vs. Conditional Multivariate Predictions for Anopheles troparvus") +
  xlab("Predictions from Univariate Probit ") + ylab("Conditional Predictions from Multivariate Probit") +
  labs( color = "PA of Culex perexiguus", shape = "PA of Anopheles troparvus") + ylim(0,1) +
  geom_abline(slope=1, intercept=0)


#We can see that the conditioning on the PA of CUlex has a clear effect on predictions of 
#mv-probit. It's roughly .125 high, if Culex is present. Prediction have roughly, but not exactly a
#slope of 1.

#####################################6.2.b Uncertainty in the predictions#######################
#We take a look at the dispersion of the predictions for the three prediction types.  
#We consider the the entire dataset to get a large sample size. No need to withold observations
#or do it on the test set, because we don't use the "true values/observations"
####For Culex perexiguus
###For the univariate model; We need to take the predictions on the probability scale
#We obtain them, getting the posterior predictions of the linear predictor (posterior_linpred)
#and transforming these value with the link function (transform = T)
pred_cpuv <- posterior_linpred(fit_cp, transform = T, newdata = df, seed = 333)

#mean of the (standard deviations of the  predictions per observation)
sd_cpuv <- apply(pred_cpuv, 2, sd) %>% mean
#sd roughly .07

#for the unconditional predictions from the multivariate model
pred_mul <- predict_mul(fit_mul, rhs_formula, df, n_max = 4000)
sd_cp_mvun <- apply(pred_mul$p1_uncond, 1, sd) %>% mean
#sd (.07) is the same as in the univariate case

#For the conditional predictions
#build a matrix of the simulations conditioned on the PA of Anopheles
pred_cpmul_con <- pred_mul$p1_cond #First take only the predictions conditioned on present Anopheles
pred_cpmul_con[df$Anatropre == 0,] <- pred_mul$p1_cond_0[df$Anatropre == 0,] #'Replace the cases in which
#Anopheles is actually absent in df with the predictions conditioned on absent Anopheles 

#calculate the mean of the sds of the predictions per observation
sd_cp_mvco <- apply(pred_cpmul_con, 1, sd) %>% mean
#sd=.07 >> same as in the other two cases

#Plotting the SDs of predictions per observation against each other (univariate vs. unconditional vs.
#vs. conditional predictions)
#'making the ggplot dataframe with all the standard deviations, the PA of Anopheles and variable
#'indicating whether observation was part of the training or test set
gg_sd_cp <- data.frame(apply(pred_cpuv, 2, sd), apply(pred_mul$p1_uncond, 1, sd),
                       apply(pred_cpmul_con, 1, sd), df$Anatropre, seq(1:nrow(df)) %in% train_id)
names(gg_sd_cp) <- c("uv", "mvun", "mvco", "at", "train")
#make factor out of at
gg_sd_cp$at <- factor(gg_sd_cp$at, levels = c(0, 1), labels = c("absent", "present"))
gg_sd_cp$train <- factor(gg_sd_cp$train, levels = c(True, False), labels = c("train", "test"))

#Univariate vs. unconditional multivariate
ggplot(gg_sd_cp, aes(x = uv, y = mvun)) + geom_point() + geom_abline(slope = 1, intercept = 0) +
  ylim(0,.2) + xlim(0, .2) + labs(title = "Standard Deviations of Univariate vs. Unconditional Multivariate Predictions \nper Observation for Culex perexiguus",
                                  y ="Standard Deviation of Multivariate Unconditional Predictions",
                                  x = "Standard Deviation of Univariate Predictions")
#points are on the identity line >> SDs of are the same across observations

#Univariate vs. conditional multivariate
ggplot(gg_sd_cp, aes(x = uv, y = mvco, color = at)) + geom_point() + geom_abline(slope = 1, intercept = 0) +
  ylim(0,.2) + xlim(0, .2) + labs(title = "Standard Deviations of Univariate and Conditional Multivariate Predictions per Observation for Culex perexiguus",
                                  y ="Standard Deviation of Multivariate Conditional Predictions",
                                  x = "Standard Deviation of Univariate Predictions", color = "PA of Anopheles atroparvus")

#Univariate vs. conditional multivariate +indication of set (train or test)
ggplot(gg_sd_cp, aes(x = uv, y = mvco, color = train)) + geom_point() + geom_abline(slope = 1, intercept = 0) +
  ylim(0,.2) + xlim(0, .2) + labs(title = "Standard Deviations of Univariate and Conditional Multivariate Predictions per Observation for Culex perexiguus",
                                  y ="Standard Deviation of Multivariate Conditional Predictions",
                                  x = "Standard Deviation of Univariate Predictions", color = "Set")

#'Interesting result: For a certain range of SDs, the multivariate model has larger SDs than the
#'univariate model and vice versa. Moreover, for more points the standard deviations of univariate
#'predictions are higher than for multivariate predictions, but if sd of multivariate is higher than
#'the difference is larger than than in the case of higher SDs for univariate predictions.

#Plotting uncertainties of predictions with boxplots for single randomly drawn observations
#draw the observation
i <- sample(seq_len(nrow(df)), size = 1)

#make the dataframe with all the predictions of observation i for all three prediction types
d_gg_bpso <- data.frame(c(pred_cpuv[,i], pred_mul$p1_uncond[i,], pred_cpmul_con[i,]))
names(d_gg_bpso) <- c("pred")
#adding a column with the prediction types
d_gg_bpso$type <- c(rep("Univariate", length(pred_cpuv[,i])),
                    rep("Unconditional Multivariate", length(pred_mul$p1_uncond[i,])),
                    rep("Conditional Multivariate", length(pred_cpmul_con[i,])))
#adding a column with the observation #
d_gg_bpso$obs <- i
# Boxplot shows Boxes: 25% and 75 % Quantile; vertical line: median;
#lower whisker = smallest observation greater than or equal to lower hinge - 1.5 * IQR; 
#points: "remaining" outliers 
ggplot(d_gg_bpso, aes(x = pred, y = type)) +
  geom_boxplot(aes(color = type)) + labs(title = paste0("Boxplot of 4000-simulated Predictions for Observation ", i),
                                         y = "Prediction Type",
                                         subtitle = "Culex perexiguus",
                                         x = "Prediction on the Probability Scale") + xlim(0,1) +
  theme(legend.position = "none") 
#'>>Univariate and unconditional are very alike. The conditional multivariate usuallaly has a
#'different median and sometimes also different width of the box.

####95%-Confidence bands of predictions along predicted probabilities of all three predictions
#types

#calculate the .025 and .975 quantiles for each prediction type and each observation on the
#simulated predictions
quant_cpuv <- colQuantiles(pred_cpuv, probs=c(.025, .975))
quant_cpun <- rowQuantiles(pred_mul$p1_uncond, probs=c(.025, .975))
quant_cpco <- rowQuantiles(pred_cpmul_con, probs=c(.025, .975))

#make the data frame for the "credibility intervall"confidence band" (correct word choice?)
#(ribbon argument of ggplot), the single argument is the median of the predictions
frame_cp <- data.frame(lower = c(quant_cpuv[,1], quant_cpun[,1], quant_cpco[,1]),
                       upper = c(quant_cpuv[,2], quant_cpun[,2], quant_cpco[,2]),
                       single = c(colMedians(pred_cpuv), rowMedians(pred_mul$p1_uncond),
                                  rowMedians(pred_cpmul_con)),
                       type = factor(c(rep("Univariate Predictions", length(pred_cpuv)),
                                       rep("Unconditional Multivariate Predictions", length(pred_cpuv)),
                                       rep("Conditional Multivariate Predictions", length(pred_cpuv))),
                                     levels = c("Univariate Predictions", "Unconditional Multivariate Predictions",
                                                "Conditional Multivariate Predictions")))

#make the plot
ggplot(frame_cp, aes(single, single))+
  geom_line() +
  geom_ribbon(data=frame_cp,aes(ymin=lower,ymax=upper),alpha=0.3) +
  facet_wrap(~type, nrow =3)+ 
  labs(title = "95%-Confidence Bands of Predictions", subtitle = "Culex perexiguus",
       y = "Prediction", x = "Prediction")

#Confidence Bands look very similiar

####For Anopheles####
#Predictions for univariate model on entire dataset
pred_atuv <- posterior_linpred(fit_at, transform = T, newdata = df, seed = 333)

#mean of the (standard deviations of the  predictions per observation)
sd_atuv <- apply(pred_atuv, 2, sd) %>% mean
#sd roughly .06

#for the unconditional predictions from the multivariate model
sd_at_mvun <- apply(pred_mul$p2_uncond, 1, sd) %>% mean
#sd (.06) is the same as in the univariate case

#For the conditional predictions
#build a matrix of the simulations conditioned on the PA of Culex
pred_atmul_con <- pred_mul$p2_cond #First take only the predictions conditioned on present Culex
pred_atmul_con[df$Cxperpre == 0,] <- pred_mul$p2_cond_0[df$Cxperpre == 0,] #'Replace the cases in which
#Culex is actually absent in df with the predictions conditioned on absent Culex 

#calculate the mean of the sds of the predictions per observation
sd_at_mvco <- apply(pred_atmul_con, 1, sd) %>% mean
#sd=.07 >> roughly the same as in the other two cases

#Plotting the SDs of predictions per observation against each other (univariate vs. unconditional vs.
#vs. conditional predictions)
#'making the ggplot dataframe with all the standard deviations, the PA of Culex and variable
#'indicating whether observation was part of the training or test set
gg_sd_at <- data.frame(apply(pred_atuv, 2, sd), apply(pred_mul$p2_uncond, 1, sd),
                       apply(pred_atmul_con, 1, sd), df$Cxperpre, seq(1:nrow(df)) %in% train_id)
names(gg_sd_at) <- c("uv", "mvun", "mvco", "cp", "train")
#make factor out of at
gg_sd_at$cp <- factor(gg_sd_at$cp, levels = c(0, 1), labels = c("absent", "present"))
gg_sd_at$train <- factor(gg_sd_at$train, levels = c(True, False), labels = c("train", "test"))

#Univariate vs. unconditional multivariate
ggplot(gg_sd_at, aes(x = uv, y = mvun)) + geom_point() + geom_abline(slope = 1, intercept = 0) +
  ylim(0,.2) + xlim(0, .2) + labs(title = "Standard Deviations of Univariate vs. Unconditional Multivariate Predictions \nper Observation for Anopheles atroparvus",
                                  y ="Standard Deviation of Multivariate Unconditional Predictions",
                                  x = "Standard Deviation of Univariate Predictions")
#points are on the identity line >> SDs of are the same across observations

#Univariate vs. conditional multivariate
ggplot(gg_sd_at, aes(x = uv, y = mvco, color = cp)) + geom_point() + geom_abline(slope = 1, intercept = 0) +
  ylim(0,.2) + xlim(0, .2) + labs(title = "Standard Deviations of Univariate and Conditional Multivariate Predictions per Observation for Anopheles troparvus",
                                  y ="Standard Deviation of Multivariate Conditional Predictions",
                                  x = "Standard Deviation of Univariate Predictions", color = "PA of for Culex perexiguus")

#Univariate vs. conditional multivariate +indication of set (train or test)
ggplot(gg_sd_at, aes(x = uv, y = mvco, color = train)) + geom_point() + geom_abline(slope = 1, intercept = 0) +
  ylim(0,.2) + xlim(0, .2) + labs(title = "Standard Deviations of Univariate and Conditional Multivariate Predictions per Observation",
                                  y ="Standard Deviation of Multivariate Conditional Predictions",
                                  x = "Standard Deviation of Univariate Predictions", color = "Set")

#'Interesting result: For a certain range of SDs, the multivariate model has larger SDs than the
#'univariate model and vice versa. Moreover, for more points the standard deviations of univariate
#'predictions are higher than for multivariate predictions, but if sd of multivariate is higher than
#'the difference is larger than than in the case of higher SDs for univariate predictions.

#Plotting uncertainties of predictions with boxplots for single randomly drawn observations
#draw the observation
i <- sample(seq_len(nrow(df)), size = 1)

#make the dataframe with all the predictions of observation i for all three prediction types
d_gg_bpso <- data.frame(c(pred_atuv[,i], pred_mul$p2_uncond[i,], pred_atmul_con[i,]))
names(d_gg_bpso) <- c("pred")
#adding a column with the prediction types
d_gg_bpso$type <- c(rep("Univariate", length(pred_atuv[,i])),
                    rep("Unconditional Multivariate", length(pred_atuv[,i])),
                    rep("Conditional Multivariate", length(pred_atuv[,i])))
#adding a column with the observation #
d_gg_bpso$obs <- i
# Boxplot shows Boxes: 25% and 75 % Quantile; vertical line: median;
#lower whisker = smallest observation greater than or equal to lower hinge - 1.5 * IQR; 
#points: "remaining" outliers 
ggplot(d_gg_bpso, aes(x = pred, y = type)) +
  geom_boxplot(aes(color = type)) + labs(title = paste0("Boxplot of 4000-simulated Predictions for Observation ", i),
                                         y = "Prediction Type",
                                         subtitle = "Culex perexiguus",
                                         x = "Prediction on the Probability Scale") + xlim(0,1) +
  theme(legend.position = "none") 
#'>>Univariate and unconditional are very alike. The conditional multivariate usuallaly has a
#'different median and sometimes also different width of the box.

####95%-Confidence bands of predictions along predicted probabilities of all three predictions
#types

#calculate the .025 and .975 quantiles for each prediction type and each observation on the
#simulated predictions
quant_atuv <- colQuantiles(pred_atuv, probs=c(.025, .975))
quant_atun <- rowQuantiles(pred_mul$p2_uncond, probs=c(.025, .975))
quant_atco <- rowQuantiles(pred_atmul_con, probs=c(.025, .975))

#make the data frame for the "credibility intervall"confidence band" (correct word choice?)
#(ribbon argument of ggplot), the single argument is the median of the predictions
frame_at <- data.frame(lower = c(quant_atuv[,1], quant_atun[,1], quant_atco[,1]),
                       upper = c(quant_atuv[,2], quant_atun[,2], quant_atco[,2]),
                       single = c(colMedians(pred_atuv), rowMedians(pred_mul$p2_uncond),
                                  rowMedians(pred_atmul_con)),
                       type = factor(c(rep("Univariate Predictions", length(pred_atuv)),
                                       rep("Unconditional Multivariate Predictions", length(pred_atuv)),
                                       rep("Conditional Multivariate Predictions", length(pred_atuv))),
                                     levels = c("Univariate Predictions", "Unconditional Multivariate Predictions",
                                                "Conditional Multivariate Predictions")))

#make the plot
ggplot(frame_at, aes(single, single))+
  geom_line() +
  geom_ribbon(data=frame_at,aes(ymin=lower,ymax=upper),alpha=0.3) +
  facet_wrap(~type, nrow =3)+ 
  labs(title = "95%-Confidence Bands of Predictions", subtitle = "Anpheles troparvus",
       y = "Prediction", x = "Prediction")

#Confidence Bands look very similiar
