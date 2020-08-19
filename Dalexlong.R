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
library(tibble)
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
#IMPORTANT: CHECK B-parameter of model_parts and adjust it! Maybe to 50... 
#I just chose 5 to fasten calculations, Moreover check ng-argument for all 
#the gjam simulations!
#Define our covariates as xdata
xdata <- train[,c("Mes", "IA_500", "NDVIBEF_2000")]
# create custom predict function for rstanarm "posterior_predict function
pred_uv <- function(model, newdata)  {
  return(posterior_predict(model, newdata) %>% apply(2, mean))
}
#makes a big difference, whether you take mean or median!, I take the mean
#bc we do it the same way in gjam
#create the explain object (core object of DALEX) which contains the data and
#the predict function
dal_cp <- explain(fit_fin_cp, xdata, y = train$Cxperpre, predict_function = pred_uv,
                  type = "classification", label = "Univariate Probit")


# calculate the permutation-based variable importance measure (the difference
#(1-AUC) between original data and permuted data per covariate)
set.seed(1980)
vi_cp <- model_parts(dal_cp, type = "difference", B = 50)
plot(vi_cp) + 
  labs(title = "Variable Importance over 50 Permuations", subtitle = "created for the univariate probit model of Culex perexiguus") 
#So, Mes is the  most important variable, followed by IA_500 and at last 
#NDVIBEF_2000. When permuting the Mes variable entries and then predicting
#our response, the resulting AUC is roughly .12 worse than for the 
#predictions without the permutations. 


#####For gjam unconditional predictions

########You need to add your quadratic terms and shit >> maybe it makes sense to do this in your
#script after you selected your final model, remember you already did this with your response 
#curves!


# create custom predict function, this is a little trickier, because we need
#to feed it xdata with the three columns (Mes, IA_500 and NDVIBEF) (bc we 
#want to get the variable importance of these three variables), but for
#gjamPredict we need to modify the data so that the function works.
pred_gj <- function(model, newdata)  {
  #prepare data for prediction
  #convert the factor mes to dummies
  newdata <- gjamDummy(newdata$Mes, newdata)
  #add quadratic terms
  newdata$"I(IA500^2)" <- (newdata$IA_500)^2
  newdata$"I(NDVIBEF2000^2)" <- (newdata$NDVIBEF_2000)^2
  
  #make a newdata list for gjam
  newdata <- list(xdata = newdata, nsim = 4000)
  #Doing the predictions
  pre <- gjamPredict(output = model, newdata = newdata)
  #return the mean of the y-chains for every observation
  return(pre$sdList$yMu[,1])
}

#make the explainer
dal_cpgj <- explain(joint_fin, xdata, y = train$Cxperpre,
                    predict_function = pred_gj, type = "classification",
                    label = "Unconditional Multivariate Probit")

#permutation-based variable importance measure
set.seed(1980)
vi_cpgj <- model_parts(dal_cpgj, type = "difference", B = 50) # du musst auf jeden Fall number of samples
#reduzieren >> dauert nämlich forever

#plot the results
plot(vi_cpgj) + labs(title = "Variable Importance over 50 Permutations", subtitle = "created for the multivariate probit model of Culex perexiguus")
#Most important variable is IA_500, then Mes and then NDVIBEF_2000

####################Doing it for conditional predictions


# create custom predict function, adding the presence-absence of Anopheles
#as ydataCond to newdata list, but treating this condition NOT as a covariate
#Probably not that important, maybe just for the "comparison plot"
pred_gjco <- function(model, newdata)  {
  #prepare data for prediction
  #convert the factor mes to dummies
  newdata <- gjamDummy(newdata$Mes, newdata)
  #add quadratic terms
  newdata$"I(IA500^2)" <- (newdata$IA_500)^2
  newdata$"I(NDVIBEF2000^2)" <- (newdata$NDVIBEF_2000)^2
  #make a newdata list for gjam including Anopheles' PA
  newdata <- list(xdata = newdata, ydataCond = y_train[,2], nsim = 4000)
  #Doing the predictions
  pre <- gjamPredict(output = model, newdata = newdata)
  #return the mean of the y-chains for every observation
  return(pre$sdList$yMu[,1])
}

#make the explainer
dal_cpgjco <- explain(joint_fin, xdata, y = train$Cxperpre,
                      predict_function = pred_gjco, type = "classification",
                      label = "Conditional Multivariate Probit")

#permutation-based variable importance measure
set.seed(1980)
#DOing it for 50 permutations
vi_cpgjco <- model_parts(dal_cpgjco, type = "difference", B = 50) # du musst auf jeden Fall number of samples

#plots the results
plot(vi_cpgjco) + labs(title = "Variable Importance", subtitle = "created for the multivariate probit model of Culex perexiguus")
#differences are much smaller, because the "conditioning" can make up for most
#of the predictive power loss due to "permuting" the specific covariate

#plot all three: 
#adding the other prediction types
d_gg <- tibble(var = c(vi_cp$variable, vi_cpgj$variable, vi_cpgjco$variable), loss = c(vi_cp$dropout_loss, vi_cpgj$dropout_loss, vi_cpgjco$dropout_loss), type = c(vi_cp$label, vi_cpgj$label, vi_cpgjco$label))
#converting var and type to a factor
d_gg$var <- as.factor(d_gg$var)
d_gg$type <- factor(d_gg$type, levels = c("Univariate Probit", "Unconditional Multivariate Probit",
                                          "Conditional Multivariate Probit"))

#deleting baseline and full model 
d_gg <-d_gg[!(d_gg$var == "_full_model_" | d_gg$var == "_baseline_"),]


ggplot(d_gg, aes(x = loss, y = var)) +
  geom_boxplot(aes(color = type)) + # Boxplot shows Boxes: 25% and 75 % Quantile; vertical line: median; lower whisker = smallest observation greater than or equal to lower hinge - 1.5 * IQR
  facet_wrap(~type, nrow =3) +
  labs(title = "Variable Importance Boxplots",
       subtitle = "Created for multi- and univariate models of Culex perexiguus with 50 permutation runs", y = "Variable",
       x = "(1 - AUC)-Loss after Permutations") +
  # Suppress the legend since color isn't actually providing any information
  theme(legend.position = "none") 


###Treat the conditioning as a covariate (inlcude it as one variable in your
#variable importance analysis)


#xdata with PA of Anopheles
xdata_con <- as_tibble(cbind(xdata, y_train[,2]))

#make the predict function, including PA-anopheles as a covariate to
#permute/manipulate
pred_con <- function(model, newdata)  {
  #prepare data for prediction
  #convert the factor mes to dummies
  newdata <- gjamDummy(newdata$Mes, newdata) %>% as_tibble()
  #add quadratic terms
  newdata$"I(IA500^2)" <- (newdata$IA_500)^2
  newdata$"I(NDVIBEF2000^2)" <- (newdata$NDVIBEF_2000)^2
  #make a newdata list for gjam
  newdata <- list(xdata = newdata, ydataCond = newdata[,match("Anatropre", names(newdata))], nsim = 4000)
  #Doing the predictions
  pre <- gjamPredict(output = model, newdata = newdata)
  #return the mean of the y-chains for every observation
  return(pre$sdList$yMu[,1])
}

#make the explainer
dal_cp_con <- explain(joint_fin, xdata_con, y = train$Cxperpre,
                      predict_function = pred_con, type = "classification")

#permutation-based variable importance measure
set.seed(1980)
#DOing it for 10 permutations
vi_cp_con <- model_parts(dal_cp_con, type = "difference", B = 50)

#plotting the results
plot(vi_cp_con) + labs(title = "Variable Importance", subtitle = "created for the multivariate probit model of Culex perexiguus conditional on Anopheles troparvus")
#Anatrope label muss noch geändert werden
#Antrope ist die wichtigste Variable...


###########################################################################
###For Anopheles
#make the explainer
dal_at <- explain(fit_fin_at, xdata, y = train$Anatropre, predict_function = pred_uv,
                  type = "classification", label = "Univariate Probit")


# calculate the permutation-based variable importance measure 
set.seed(1980)
vi_at <- model_parts(dal_at, type = "difference", B = 50)
#plot it
plot(vi_at) + 
  labs(title = "Variable Importance over 50 Permuations", subtitle = "created for the univariate probit model of Anopheles troparvus") 
#So, Mes is the  most important variable, followed by NDVIBEF_2000 and at last 
#IA_500. When permuting the Mes variable entries and then predicting
#our response, the resulting AUC is roughly .18 worse than for the 
#predictions without the permutations. 

######For unconditional predictions in gjam

#We need to slightly change our predict function (just the last line)
pred_gjat <- function(model, newdata)  {
  #prepare data for prediction
  #convert the factor mes to dummies
  newdata <- gjamDummy(newdata$Mes, newdata)
  #add quadratic terms
  newdata$"I(IA500^2)" <- (newdata$IA_500)^2
  newdata$"I(NDVIBEF2000^2)" <- (newdata$NDVIBEF_2000)^2
  
  #make a newdata list for gjam
  newdata <- list(xdata = newdata, nsim = 4000)
  #Doing the predictions
  pre <- gjamPredict(output = model, newdata = newdata)
  #return the mean of the y-chains for every observation, this time of at
  return(pre$sdList$yMu[,2])
}

#make the explainer
dal_atgj <- explain(joint_fin, xdata, y = train$Anatropre,
                    predict_function = pred_gjat, type = "classification", 
                    label = "Unconditional Multivariate Probit")

#permutation-based variable importance measure
set.seed(1980)
vi_atgj <- model_parts(dal_atgj, type = "difference", B = 50) 
#pretty similar to univariate model!

#plot the results
plot(vi_atgj) + labs(title = "Variable Importance over 50 Permutations", subtitle = "created for the multivariate probit model of Anopheles troparvus")
#Most important variable is Mes, then NDVI_BEF and then IA, similar results

####################Doing it for conditional predictions


# create custom predict function, adding the presence-absence of Anopheles
#as ydataCond to newdata list, but treating this condition NOT as a covariate
#Probably not that important, maybe just for the "comparison plot"
pred_atgjco <- function(model, newdata)  {
  #prepare data for prediction
  #convert the factor mes to dummies
  newdata <- gjamDummy(newdata$Mes, newdata)
  #add quadratic terms
  newdata$"I(IA500^2)" <- (newdata$IA_500)^2
  newdata$"I(NDVIBEF2000^2)" <- (newdata$NDVIBEF_2000)^2
  #make a newdata list for gjam including Culex' PA
  newdata <- list(xdata = newdata, ydataCond = y_train[,1], nsim = 4000)
  #Doing the predictions
  pre <- gjamPredict(output = model, newdata = newdata)
  #return the mean of the y-chains for every observation
  return(pre$sdList$yMu[,2])
}

#make the explainer
dal_atgjco <- explain(joint_fin, xdata, y = train$Anatropre,
                      predict_function = pred_atgjco, type = "classification",
                      label = "Conditional Multivariate Probit")

#permutation-based variable importance measure
set.seed(1980)
#DOing it for 50 permutations
vi_atgjco <- model_parts(dal_atgjco, type = "difference", B = 50) 

#plots the results
plot(vi_atgjco) + labs(title = "Variable Importance", subtitle = "created for the multivariate probit model of Anopheles troparvus")
#differences are much smaller, because the "conditioning" can make up for most
#of the predictive power loss due to "permuting" the specific covariate

#######plot all three: 
#adding the other prediction types
d_ggat <- tibble(var = c(vi_at$variable, vi_atgj$variable, vi_atgjco$variable),
                 loss = c(vi_at$dropout_loss, vi_atgj$dropout_loss, vi_atgjco$dropout_loss),
                 type = c(vi_at$label, vi_atgj$label, vi_atgjco$label))
#converting var and type to a factor
d_ggat$var <- as.factor(d_ggat$var)
d_ggat$type <- factor(d_ggat$type, levels = c("Univariate Probit", "Unconditional Multivariate Probit",
                                          "Conditional Multivariate Probit"))

#deleting baseline and full model 
d_ggat <-d_ggat[!(d_ggat$var == "_full_model_" | d_ggat$var == "_baseline_"),]


ggplot(d_ggat, aes(x = loss, y = var)) +
  geom_boxplot(aes(color = type)) + # Boxplot shows Boxes: 25% and 75 % Quantile; vertical line: median; lower whisker = smallest observation greater than or equal to lower hinge - 1.5 * IQR
  facet_wrap(~type, nrow =3) +
  labs(title = "Variable Importance Boxplots",
       subtitle = "Created for multi- and univariate models of Anopheles troparvus with 50 permutation runs", y = "Variable",
       x = "(1 - AUC)-Loss after Permutations") +
  # Suppress the legend since color isn't actually providing any information
  theme(legend.position = "none") 


###Treat the conditioning as a covariate (inlcude it as one variable in your
#variable importance analysis)

#xdata with PA of Culex
xdata_con <- as_tibble(cbind(xdata, y_train[,1]))

#make the predict function, including PA-Culex as a covariate to
#permute/manipulate
pred_atcon <- function(model, newdata)  {
  #prepare data for prediction
  #convert the factor mes to dummies
  newdata <- gjamDummy(newdata$Mes, newdata) %>% as_tibble()
  #add quadratic terms
  newdata$"I(IA500^2)" <- (newdata$IA_500)^2
  newdata$"I(NDVIBEF2000^2)" <- (newdata$NDVIBEF_2000)^2
  #make a newdata list for gjam
  newdata <- list(xdata = newdata, ydataCond = newdata[,match("Cxperpre", names(newdata))], nsim = 4000)
  #Doing the predictions
  pre <- gjamPredict(output = model, newdata = newdata)
  #return the mean of the y-chains for every observation
  return(pre$sdList$yMu[,2])
}

#make the explainer
dal_at_con <- explain(joint_fin, xdata_con, y = train$Anatropre,
                      predict_function = pred_atcon, type = "classification")

#permutation-based variable importance measure
set.seed(1980)
#DOing it for 10 permutations
vi_at_con <- model_parts(dal_at_con, type = "difference", B = 50)

#plotting the results
plot(vi_at_con) + labs(title = "Variable Importance", subtitle = "created for the multivariate probit model of Anopheles troparvus conditional on Culex perexiguus")
#Anatrope label muss noch geändert werden
#Culex isnt that important o_o; is that in accordance with the AUCs?

###########################################################################
#Check mal AUC von beiden complete Models!
#mach vielleicht auch beides in einen plot, dann sieht man direkt, wie die
#Modelle zu einander stehen
auc_cpgj <- loss_one_minus_auc(observed = train$Cxperpre,
                               predicted = pred_gj(joint_fin, xdata))
auc_cpuv <- loss_one_minus_auc(observed = train$Cxperpre,
                               predicted = pred_uv(fit_fin_cp, xdata))
auc_atuv <- loss_one_minus_auc(observed = train$Anatropre,
                               predicted = pred_uv(fit_fin_at, xdata))
#>> are very similiar :), like we hypothesized

