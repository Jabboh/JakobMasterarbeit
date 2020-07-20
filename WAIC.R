#stan_glm function macht Faxen, der "initialization  failed"-error macht die Probleme

rm(list=ls())
setwd("C:\\Users\\jakob\\Documents\\JakobMasterarbeit\\Data")
#loading packages
library(readxl) #to read in the data
library(rstanarm) #Doing Bayesian probit GLMs for single species
library(pROC) #calculating the AUC
library(gjam) #Doing the joint estimation with the GJAM modelling aproach
library(ggplot2) #for plotting
####Data Preperation

#read in the data (Monthly species PA data for seven different mosquito species
#and according environmental covariates)
df <- read_excel("MonthlyData.xlsx")
str(df)
#We ignore An_atroparvus, because it does not seem to be PA-data.
#It looks like abundance data

#extract the PA data for all species
spec <- df[,7:14]
#deleting the An_atroparvus column
spec[,"An_atroparvus"] <- NULL
#To select two species for a first analysis, we look for species that roughly occur in 50%
#of the observations.
summary(spec)
#Hence, we select Cxperpre and Anatropre for our first analysis:
y <- spec[,c("Cxperpre", "Anatropre")]
#For these two species, Roiz et al. use the following covariates to explain their 
#monthly PAs: Inundation area (500 m buffer), NDVI (500 m buffer), Inundation area (2000 m)
#and NDVI month before (2000 m buffer). Since NDVI 500 m buffer and NDVI 2000 m buffer
#will be highly correlated, I will leave out NDVI 2000 m buffer.

#Normalizing Covariates: Otherwise, interaction terms hardy interpretable and skewed
df[,17:ncol(df)] <- scale(df[,17:ncol(df)])

#Split the data set in training (70 %) and test (30%) set
train_size <- floor(0.7 * nrow(y))
#set a random seed for replicability
set.seed(333)
#sample the training IDs
train_id <- sample(seq_len(nrow(y)), size = train_size)

#partition data into train and test set
train <- df[train_id, ]
test <- df[-train_id, ]

#####Doing the probit regressions in a Bayesian setting for every species seperately
#For this we use the stan_glm function from the rstanarm-package

#fitting the model for Culex Perexiguus
fit_cp <- stan_glm(Cxperpre ~ (IA_500 + NDVI_500 + NDVIBEF_2000)^2 + I(IA_500^2) + I(NDVI_500^2) + I(NDVIBEF_2000^2), data = train, family = binomial(link = "probit"), init_r = 1, seed = 333)
summary(fit_cp)

#WAIC: http://mc-stan.org/rstanarm/reference/loo.stanreg.html

waic(fit_cp)
loo(fit_cp)


####Fitting a joint model of both species with GJAM

#Define the model settings
types <- c("PA", "PA") #variable types of the responses
s <- length(types) # number of species
#y-data to train the model
y_train <- y[train_id, ]

#define model/algorithm parameters
ml   <- list(ng = 2000, burnin = 100, typeNames = types)
#is 2000 for number of Gibbs steps ok?

#runnig GJAM
joint <- gjam(~ (IA_500 + NDVI_500 + NDVIBEF_2000)^2 + I(IA_500^2) + I(NDVI_500^2) + I(NDVIBEF_2000^2), ydata = y_train, xdata = train, modelList = ml)
summary(joint)
joint$chains
