rm(list=ls())

#loading packages
library(readxl) #to read in the data
library(rstanarm) #Doing Bayesian probit GLMs for single species
library(pROC) #calculating the AUC
library(gjam) #Doing the joint estimation with the GJAM modelling aproach

#read in the data (Species PA Data and covariates)
df <- read_excel("MonthlyData.xlsx")
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
library(jtools) #helpful plotting tools to compare coefficients
#https://cran.r-project.org/web/packages/jtools/vignettes/summ.html
###Conditional Prediction
#Culex perexiguus conditioned on Anopheles atroparvus

#Anopheles atroparvus conditioned on Culex perexiguus