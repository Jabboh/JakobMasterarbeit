#############################Results of my analyis
#1. Data preparation and fitting of final model
#6. Results
#6.1. In-sample (not on test data):
#a. Comparing the coefficients: Size and credibility intervals
#b. Correlations between the responses / Residual Correlation Parameter of gjam
#c. Response Curves ("in-s
#d. Variable importance
#6.2. Out-of-Sample
#a. Conditional Predictions of gjam vs. Unconditional Predictions of 
#gjam vs. predictions of univariate model
#b. Comparing the uncertainty of the different "prediction"-types
rm(list=ls())
setwd("C:\\Users\\jakob\\Documents\\JakobMasterarbeit\\Data")
#install.packages("gridExtra")

#loading packages
library(readxl) #to read in the data
library(rstanarm) #Doing Bayesian probit GLMs for single species
library(pROC) #calculating the AUC
library(gjam) #Doing the joint estimation with the GJAM modelling aproach
library(ggplot2) #for plotting
library(DHARMa) # for checking in-sample validity of the models
library(dplyr) # for simplified syntax and neater code
library(corrplot)
library(loo) #to calculate WAIC
library(bayesplot) #Some handy features for plotting in the Bayesian realm
library(gridExtra) # for plotting multiple graphs in one window
library(DALEX) #Explanatory Model Analysis / Variable Importance
#### 1.Data Preperation
#read in the data (Monthly species PA data for seven different mosquito species
#and according environmental covariates)
df <- read_excel("MonthlyData.xlsx")

#Checking the rough structure of the data set
str(df)
summary(df$Fecha) # Dates in 2006 dont make any sense. I assume that they put by
#accident 2006 instead of 2010. So I change these dates
df$Fecha <- as.POSIXct(sub("2006", "2010", df$Fecha))

#Transform "Mes" (month when record was taken) into a factor variable
df$Mes <- factor(df$Mes, levels =c("Abril", "Mayo", "Junio", "Julio", "Agosto",
                                   "Septiembre"))

#Reading in the spatial coordinates of the different trap locations
coords <- read_excel("Traps_coordenadas_geograficas.xls")

#Trap and Area mean the same thing, so we change the name in df from "area" to "trap"
names(df)[names(df)=="Area"] <- "trap"
#Canada is spelt differently, so I change the spelling in df to "Cañada"
df[,"trap"] <- lapply(df[,"trap"], gsub, pattern = "Ca?da", replacement = "Cañada",
                      fixed = T)

#adding lon-lat column to the data frame df
df <- merge(df, coords[, c("trap", "Norte", "Oeste")], by="trap", all.x= T, sort = F)

#Selecting our two species (Perexiguus & Anopheles troparvus) for the analysis. Our
#rationale is to select the two species that roughly occur in 50 % of the observations
#in order to get the maximum variation. We ignore the column "An_atroparvus", 
#because it does not seem to be PA-data (rather abundance data).

#extract the PA data for all species
spec <- df[,7:14]
#deleting the An_atroparvus column
spec[,"An_atroparvus"] <- NULL
#Taking a look at occurence rate across all observations (the Mean)
summary(spec)
#Hence, we select Cxperpre and Anatropre for our first analysis:
y <- spec[,c("Cxperpre", "Anatropre")]

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

#Fitting final univariate models

#for Culex perexiguus

fit_fin_cp <- stan_glm(Cxperpre ~ Mes + IA_500 + NDVIBEF_2000 + I(IA_500^2) + I(NDVIBEF_2000^2),
                       data = train, family = binomial(link = "probit"), init_r = 1.4,  seed = 333)

#for anopheles
fit_fin_at <- stan_glm(Anatropre ~ Mes + IA_500 + NDVIBEF_2000 + I(IA_500^2) +  I(NDVIBEF_2000^2),
                       data = train, family = binomial(link = "probit"), init_r = 1.4, seed = 333)


#c. Preparing Data for gjam

#Preparing the month factor for gjam (creating #factor-dummy variables) with our gjamDummy function
#first argument: vector of the factor variable; second argument: dataframe to which new dummies
#should be added
source("gjamDummy.R")
data_dum <- gjamDummy(df$Mes, df)
#subset training data set
train_gj <- data_dum [train_id, ]
#Define the model settings
types <- c("PA", "PA") #variable types of the responses
s <- length(types) # number of species
ml   <- list(ng = 4000, burnin = 1000, typeNames = types)
#runnig GJAM
joint_fin <- gjam(~  Mayo + Junio + Julio + Agosto + Septiembre + IA_500 +
                    NDVIBEF_2000 + I(IA_500^2) + I(NDVIBEF_2000^2),
                  ydata = y_train, xdata = train_gj, modelList = ml)

##############################Dalex for univariate model
#xdata
xdata <- train[,c("Mes", "IA_500", "NDVIBEF_2000")]
# create custom predict function
pred <- function(model, newdata)  {
  return(posterior_predict(model, newdata) %>% apply(2, median))
}

dal_cp <- explain(fit_fin_cp, xdata, y = train$Cxperpre, predict_function = pred,
                  type = "classification")

#permutation-based variable importance measure
set.seed(1980)
vi_cp <- model_parts(dal_cp)
plot(vi_cp) #why label = lm??? I specified classification above?
#I, also dont understand "1 - AUC loss" o_o 

#####For gjam unconditional predictions

########You need to add your quadratic terms and shit >> maybe it makes sense to do this in your
#script after you selected your final model, remember you already did this with your response 
#curves!
xdata_gj <- train_gj[,c("IA_500", "NDVIBEF_2000", "Abril", "Mayo", "Junio", "Julio", "Agosto",
                        "Septiembre")]
xdata_gj$intercept <- rep(1 , nrow(xdata_gj))

# create custom predict function
pred <- function(model, newdata)  {
  #prepare data for predictions
  newdata <- list(xdata = newdata, nsim = 4000)
  #the predictions
  pred <- gjamPredict(output = model, newdata = newdata)
  
  return()
}