rm(list=ls())
setwd("C:\\Users\\jakob\\Documents\\JakobMasterarbeit")
#loading packages
#install.packages("bayesplot")
library(readxl) #to read in the data
library(rstanarm) #Doing Bayesian probit GLMs for single species
library(pROC) #calculating the AUC
library(gjam) #Doing the joint estimation with the GJAM modelling aproach

#read in the data (Species PA Data and covariates)
df <- read_excel("MonthlyPAData.xlsx")
str(df)
#We ignore An_atroparvus, because it does not seem to be PA-data. No clue why?
#Looks like abundance data

#extract the PA data for all species
spec <- df[,7:14]
#deleting the An_atroparvus column
spec[,"An_atroparvus"] <- NULL
#To select two species for a first analysis, we look for species that roughly occur at 50%
#of the observations
summary(spec)
#Hence, we select Cxperpre and Anatropre for our first analysis:
y <- spec[,c("Cxperpre", "Anatropre")]
#For these two species, Roiz et al. use the following covariates to explain their 
#monthly PAs: Inundation area (500 m buffer), NDVI (500 m buffer), Inundation area (2000 m)
#and NDVI month before (2000 m buffer). Since NDVI 500 m buffer and NDVI 2000 m buffer
#will be highly correlated, I will leave out NDVI 2000 m buffer.

###Split the data set in training (70 %) and test (30%) set
train_size <- floor(0.7 * nrow(y))
#set a random seed for replicability
set.seed(333)
#sample the training IDs
train_id <- sample(seq_len(nrow(y)), size = train_size)

#partition data into train and test set
train <- df[train_id, ]
test <- df[-train_id, ]
#####Doing the probit regressions in a Bayesian setting for every species seperately
#fitting the model for Culex Perexiguus
fit_cxper <- stan_glm(Cxperpre ~ IA_500 + NDVI_500 + NDVIBEF_2000, data = train, family = binomial(link = "probit"), seed = 333)
summary(fit_cxper)

#fitting the model for Anopheles atroparvus
fit_anatr <- stan_glm(Anatropre ~ IA_500 + NDVI_500 + NDVIBEF_2000, data = train, family = binomial(link = "probit"), seed = 333)
summary(fit_anatr)

#predictions on the test set
#for Culex perexiguus
pred_cxper_sin <- posterior_predict(fit_cxper, newdata = test, seed = 333)
dim(pred_cxper_sin)
#Is it correct that for each test data point 4000 draws from the posterior predictive 
#distribution are made? So I could take the average of these draws as an estimation of the 
#of the "expected" predicted y-value?
summary(pred_cxper_sin)
pred_exp_cxper_sin <- colMeans((pred_cxper_sin))

#for Anopheles atroparvus
pred_anatr_sin <- posterior_predict(fit_anatr, newdata = test, seed = 333)
pred_exp_anatr_sin <- colMeans((pred_anatr_sin))
#Quick Model Validation: To check whether model makes some sense
perf_cxper_sin <- auc(response = test$Cxperpre, predictor = pred_exp_cxper_sin)
#AUC = .74 >> seems ok for our purposes (remember: our goal is not to find the perfect
#model for our data, but rather evaluate whether knowing one species helps our predictions
#of the other one)
perf_anatr_sin <- auc(response = test$Anatropre, predictor = pred_exp_anatr_sin)
#AUC = .78 >> seems OK for our purposes

###############################Fitting a joint model of both species with GJAM and analyzing conditional (on the other species) predictions
types <- c("PA", "PA")
s <- length(types) # number of species
#y-data to ?train the model
y_train <- y[train_id, ]
#y-data to test the model
y_test <- y[-train_id, ]
#define model/algorithm parameters
ml   <- list(ng = 1000, burnin = 100, typeNames = types)
#run the model for cas and the as independent variables
joint <- gjam(~ IA_500 + NDVI_500 + NDVIBEF_2000, ydata = y_train, xdata = train, modelList = ml)
summary(joint)
###Comparison of the coefficients
#Make a table (dataframe) to store different estimates and SEs according to the models
###For Culex Perexiguus
#For the coefficients
cof_sum_px <- data.frame(matrix(ncol = 2, nrow = 4))
co <- c("Coefficients_sin", "Coefficients_gjam")
ro <- c("Intercept_sin", "IA500_sin", "NDVI500_sin", "NDVIBEF2000_sin")
colnames(cof_sum_px) <- co
rownames(cof_sum_px) <- ro
#For the SEs
se_sum_px <- data.frame(matrix(ncol = 2, nrow = 4))
co <- c("SE_sin", "SE_gjam")
colnames(se_sum_px) <- co
rownames(se_sum_px) <- ro

###For Anopheles troparvus
#For the coefficients
cof_sum_at <- cof_sum_px
se_sum_at <- se_sum_px

#Filling the tables accordingly
cof_sum_px$Coefficients_sin <-  fit_cxper$coefficients
cof_sum_px$Coefficients_gjam <- joint$parameters$betaMu[,"Cxperpre"]
cof_sum_at$Coefficients_sin <-  fit_anatr$coefficients
cof_sum_at$Coefficients_gjam <- joint$parameters$betaMu[,"Anatropre"]

se_sum_px$SE_sin <-  fit_cxper$ses
se_sum_px$SE_gjam <- joint$parameters$betaSe[,"Cxperpre"]
se_sum_at$SE_sin <-  fit_anatr$ses
se_sum_at$SE_gjam <- joint$parameters$betaSe[,"Anatropre"]

#Coefficients and SEs for Culex perexiguus
cof_sum_px
se_sum_px
#Everything looks pretty similar, as we expected (We expexted the environmental coefficients to
#be the same. Difference between coefficients way smaller than according SEs.

#Coefficients and SEs for Anopheles troparvus
cof_sum_at
se_sum_at
#Everything looks pretty similar as expected.Difference between coefficients way smaller than according SEs.


library(bayesplot) #helpful plotting tools to compare coefficients
library(ggplot2)
posterior <- as.matrix(fit_cxper)

plot_title <- ggtitle("Posterior distributions",
                      "with medians and 80% intervals")
mcmc_areas(posterior,
           pars = c("(Intercept)", "IA_500", "NDVI_500", "NDVIBEF_2000"),
           prob = 0.8) + mcmc_areas(posterior,
                                    pars = c("(Intercept)", "IA_500", "NDVI_500", "NDVIBEF_2000"),
                                    prob = 0.8) + plot_title

fit_cxper$ses

###Conditional Prediction
#Culex perexiguus conditioned on Anopheles atroparvus


#############Das hat fuktioniert

# preparing testing data set for the gjam_predict function 
test_gj <- df[-train_id,]
test_gj$intercept <- rep(1 , nrow(test_gj)) #adding the intercept
test_gj <- test_gj[,c("intercept", "IA_500", "NDVI_500", "NDVIBEF_2000")] # getting rid of all the unused variables
yy <- y[-train_id,] 

#storing input data in newdata
newdata <- list(xdata = test_gj, ydataCond = yy[,2], nsim = 200) # conditionally predict out-of-sample
#Doing the actual prediction
p2      <- gjamPredict(output = joint, newdata = newdata)
#########################

