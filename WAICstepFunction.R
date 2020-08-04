################################Analysis of Culex Perexiguus & Anopheles troparvus:
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
#6.1. In-sample:
#a. Comparing the coefficients: Size and credibility intervals
#b. Response Curves
#c. Variable importance
#d. Residual Correlation Parameter of gjam
#6.2. Out-of-Sample
#a. Conditional Predictions of gjam vs. Unconditional Predictions of 
#gjam vs. predictions of univariate model
#b. Comparing the uncertainty of the different "prediction"-types
rm(list=ls())
setwd("C:\\Users\\Jakob\\Documents\\Uni\\GCE\\Thesis\\JakobMasterarbeit\\Data")
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

####Fitting the most complex models
#The covariates Inundation Area and NDVI are available on different spatial and
#temporal scales. For the temporal scales, I include both scales (current +
#month before) in the analyis. For the spatial scale, I select one of the five
#alternatives (100 m buffer, 250 m, 500 m, 1000 m and 2000m). I decided which scale
#to choose on the basis of the models in Roiz's paper. He already has "the most
#appropriate" specifications for his Generalized Linear Mixed Models (GLMM) with a 
#binomial error distribution and binomial link function for Culex Perexiguus & 
#Anopheles troparvus. So, for example since he used "NDVI_500" for his Culex perexiguus
#model, I use also the "NDVI_500" in all my models. Justifications for my choice of 
#spatial scales:
#IA_500: In Roiz's Culex Perexiguus model
#NDVI_500: In Roiz's Culex Perexiguus model
#IABEF_2000: In contrast to IA_500, because this way we reduce the expected 
#collinearity between  IA_500 and IABEF_2000 >> "maximum variation
#NDVIBEF_2000: In Roiz's Anopheles model

#a. fitting the model for Culex Perexiguus
fit_cp <- stan_glm(Cxperpre ~ (IA_500 + NDVI_500 + IABEF_2000 + NDVIBEF_2000)^2 +
                     Mes + I(IA_500^2) + I(NDVI_500^2) + I(IABEF_2000^2) +
                     I(NDVIBEF_2000^2), data = train, refresh = 0,
                   family = binomial(link = "probit"),init_r = .7, seed = 333)

#b. fitting the model for Anopheles
fit_at <- stan_glm(Anatropre ~ (IA_500 + NDVI_500 + IABEF_2000 + NDVIBEF_2000)^2 + 
                     Mes + I(IA_500^2) + I(NDVI_500^2) + I(IABEF_2000^2) +
                     I(NDVIBEF_2000^2), data = train, refresh = 0,
                   family = binomial(link = "probit"),init_r = .7, seed = 333)

#1. Step: Removing entire covariates and their associated terms (interactions and 
#quadratic)
x <- fit_cp
#x: fitted model for which coefficients should be dropped
#y: names of the covariates ()
z <- c(names(fit_cp$coefficients)[2:5], "Mes")
step_waic_entCov <- function(x, y){
  #define the base model to which the reduced models will be compared
  base <- x
  #Create variable w_com that serves as the condition for the while-loop
  #After the first run, it is equal to the "best" model
  w_com <- matrix("start")
  rownames(w_com) <- "start"
#Create a while loop, that keeps dropping variables until, the complexer model
  #performs better than the reduced models
  while(rownames(w_com)[length(rownames(w_com))] != "base"){
    #Create an empty list to store the objects of the different models
    st <- vector("list", length(y) + 1)
    #an empty list for all the WAICs
    waics <- st
    #Run the models (removing one covariate and all its associations in each run)
    for (i in y){
      #index
      j <-match(i,y)
      #We run the model without the terms associated with y (linear, quadratic
      #and interaction terms and store the fitted model in st)
      st[[j]] <- update(base, formula=drop.terms(base$terms, grep(i, attr(x$terms, "term.labels") ), keep.response=TRUE))
      print(st[[j]]$terms)
      #Calculating the WAIC of the model and storing it in waics
      waics[[j]] <- waic(st[[j]])
    }

    #add the baseline WAIC to which the more parsimonious models will 
    #be compared (It's always the "Most" complex model at the specific point of the 
    #model-selection procedure)
    waics[[length(waics)]] <- waic(base)
    
    #a vector with all the model names (the ones without one covariate + the base)
    n <- c(y, "base")
    #naming the two storage lists accordingly
    names(st) <- n
    names(waics) <- n
    #Compare the different WAICs
    w_com <- loo_compare(waics)
    print(w_com)
    #Printing a message to the screen indicating that variable xxx will be thrown out
    print(paste0 ("The variable ", rownames(w_com)[length(n)], " and its associations will be kicked out of the model"))
    #deleting the variable that we throw out according to WAIC in the variable
    #names vector y
    y <- y[y != rownames(w_com)[length(n)]]
    #storing the model without the removed variable as the new base model
    base <- st[[which(names(st) == rownames(w_com)[length(n)])]]
  }
  #Return the final model
  return (base)#for some reason this does not work
}

####just for testing reasons (delete later)
hey <- step_waic_entCov(x,z)
st <- vector("list", length(y))
names(st) <- y
f <- list(waic(fit_cp1_ia), waic(fit_cp1_ndvi), waic(fit_cp1_mes))
apply(st, FUN = waic) 
waic(st[[2]])
