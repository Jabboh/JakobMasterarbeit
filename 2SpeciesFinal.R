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

#Checking the correlations between the covariates
x <- cbind(train$IA_500, train$NDVI_500, train$IABEF_2000, train$NDVIBEF_2000, train$Mes)
corr <- cor(x)
cor(x)
#The correlation between NDVI_500 and NDVI_2000BEF might be a problem >> Björn?
corrplot(corr, type = "upper", tl.col = "black", tl.srt = 45)
#Doing the probit regressions in a Bayesian setting for every species seperately
#For this we use the stan_glm function from the rstanarm-package

#a. fitting the model for Culex Perexiguus
fit_cp <- stan_glm(Cxperpre ~ (IA_500 + NDVI_500 + IABEF_2000 + NDVIBEF_2000)^2 +
                    Mes + I(IA_500^2) + I(NDVI_500^2) + I(IABEF_2000^2) +
                    I(NDVIBEF_2000^2), data = train,
                   family = binomial(link = "probit"),init_r = .7, seed = 333)

#b. fitting the model for Anopheles
fit_at <- stan_glm(Anatropre ~ (IA_500 + NDVI_500 + IABEF_2000 + NDVIBEF_2000)^2 + 
                    Mes + I(IA_500^2) + I(NDVI_500^2) + I(IABEF_2000^2) +
                     I(NDVIBEF_2000^2), data = train,
                   family = binomial(link = "probit"),init_r = .7, seed = 333)
#c. Fitting the multivariate gjam-model

#Preparing the month factor for gjam (creating #factor-dummy variables
#create a storage matrix
st <- data.frame(matrix(, nrow=nrow(df), ncol = 6))
#assign the factor names to the column names
names(st) <- levels(df$Mes)
#loop over all the levels
for (i in levels(df$Mes)){
  #loop over all the observations
  for (j in 1:nrow(df)){
    st[j, i] <- ifelse(df[j, "Mes"] == i, 1, 0) #fill the dummy data frame
  }
}

#combinig our original data with dummy data frame
data_dum <- cbind(df, st)
train_gj <- data_dum [train_id, ]
#Define the model settings
types <- c("PA", "PA") #variable types of the responses
s <- length(types) # number of species

#define model/algorithm parameters: 4000 gibbs steps (equivalent to the MCMCs draws
#in the rstanarm model fitting) + burnin of 1000
ml   <- list(ng = 4000, burnin = 1000, typeNames = types)

#runnig GJAM
joint <- gjam(~ (IA_500 + NDVI_500 + IABEF_2000 + NDVIBEF_2000)^2 + Mayo + Junio +
                Julio + Agosto + Septiembre + I(IA_500^2) + I(NDVI_500^2) +
                I(IABEF_2000^2) + I(NDVIBEF_2000^2), ydata = y_train, xdata = train_gj, modelList = ml)

####Internal Validation of the most complex model with DHARMa
#a. Culex Perexiguus
#Create a dharma object. For this, I specify the following arguments:
#1. 4000 simulations of fitted responses per observation
#(simulated response), 2. the observed responses (observedResponse),
# 3. the median of the 4000 simulated fitted responses which is the "expected" 
#value of the predicted y of an observation (fittedPredictedResponse),
#and other arguments concerning the specific handling of the "scaled residuals".
dharm_cp <- createDHARMa(simulatedResponse = t(posterior_predict(fit_cp)), observedResponse = fit_cp$y,
                         fittedPredictedResponse = posterior_predict(fit_cp) %>% apply(2, median), integerResponse = T, seed = 333,
                         method = "PIT")
plot(dharm_cp)
#1. QQ-plot looks good >> scaled residuals follow the expected uniform distribution
#2. KS-test also supports our assumption that the scaled residuals are uniformly 
#distributed >> our model is correctly specified.
#3. Boxplot: For both model predictions (absent or present), the scaled residuals are
#approximately equally distributed >> no bias

#Plotting the residuals against all covariates to check whether we specified the
# functional relationships correctly.

#Plot residuals against all covariates
plotResiduals(dharm_cp, train$IA_500) #looks ok, no significant problems
plotResiduals(dharm_cp, train$NDVI_500)#no significant problems 
plotResiduals(dharm_cp, train$IABEF_2000)#significant problems (.25 and .75 quantiles)
plotResiduals(dharm_cp, train$NDVIBEF_2000)#no problems
plotResiduals(dharm_cp, train$Mes) #looks good
#the quadratic forms
plotResiduals(dharm_cp, (train$IA_500)^2) # looks ok
plotResiduals(dharm_cp, (train$NDVI_500)^2) # looks ok
plotResiduals(dharm_cp, (train$IABEF_2000)^2) # looks ok
plotResiduals(dharm_cp, (train$NDVIBEF_2000)^2) # looks ok
#the interactions
plotResiduals(dharm_cp, train$IA_500*train$NDVI_500)#looks ok
plotResiduals(dharm_cp, train$IA_500*train$IABEF_2000) # looks ok
plotResiduals(dharm_cp, train$IA_500*train$NDVIBEF_2000) # looks ok
plotResiduals(dharm_cp, train$NDVI_500*train$IABEF_2000) # looks ok
plotResiduals(dharm_cp, train$NDVI_500*train$NDVIBEF_2000) # .25 is significant
plotResiduals(dharm_cp, train$IABEF_2000*train$NDVIBEF_2000)# no problems
#Overall, only two plots had significant problems. Due to the high number of plots,
#this is somewhat  expected (we performed 15*3 = 45 significance tests >> 
#hence we would expect roughly 2 significant plots (95%-level))

hist(dharm_cp)
#looks pretty flat >> thumbs up!

##Test for temporal autocorrelation
dharm_cp_auto = recalculateResiduals(dharm_cp, group = train$Fecha)
testTemporalAutocorrelation(dharm_cp_auto, time =  unique(train$Fecha))
#doesn't  seem to be a problem

#Spatial Autocorrelation

#I do not have the coordinates for trap "M29" --> one NA in the coordinates >> We need to remove that row
dharm_cp_spatial <- recalculateResiduals(dharm_cp, group = train$Norte)
testSpatialAutocorrelation(dharm_cp_spatial, 
                           x =  aggregate(train$Oeste, list(train$trap), mean)$x, 
                           y = aggregate(train$Norte, list(train$trap), mean)$x)
#No spatial autocorrelation >> Yeah!

#b. Anopheles troparvus
dharm_at <- createDHARMa(simulatedResponse = t(posterior_predict(fit_at)),
                         observedResponse = fit_at$y,
                         fittedPredictedResponse = posterior_predict(fit_at) %>% apply(2, median),
                         integerResponse = T, seed = 333, method = "PIT")
plot(dharm_at)
#1. QQ-plot looks good >> scaled residuals follow the expected uniform distribution
#2. KS-test also supports our assumption that the scaled residuals are uniformly 
#distributed >> our model is correctly specified.
#3. Boxplot: For both model predictions (absent or present), the scaled residuals are
#approximately equally distributed >> no bias

#Plotting the residuals against all covariates to check whether we specified the
# functional relationships correctly.

#Plot residuals against all covariates
plotResiduals(dharm_at, train$IA_500) #.25 quantile is significant
plotResiduals(dharm_at, train$NDVI_500)#no significant problems 
plotResiduals(dharm_at, train$IABEF_2000)#significant problem(.25 quantile)
plotResiduals(dharm_at, train$NDVIBEF_2000)#.25 quantile is significant
plotResiduals(dharm_at, train$Mes) #looks good
#the quadratic forms
plotResiduals(dharm_at, (train$IA_500)^2) # .25 quantile is significant
plotResiduals(dharm_at, (train$NDVI_500)^2) # .25 quantile is significant
plotResiduals(dharm_at, (train$IABEF_2000)^2) # .25 quantile is significant
plotResiduals(dharm_at, (train$NDVIBEF_2000)^2) # .25 quantile is significant
#the interactions
plotResiduals(dharm_at, train$IA_500*train$NDVI_500)#.25 quantile is significant
plotResiduals(dharm_at, train$IA_500*train$IABEF_2000) # .25 quantile is significant
plotResiduals(dharm_at, train$IA_500*train$NDVIBEF_2000) # .25 quantile is significant
plotResiduals(dharm_at, train$NDVI_500*train$IABEF_2000) # .25 and .5 quantile is significant
plotResiduals(dharm_at, train$NDVI_500*train$NDVIBEF_2000) # .25 is significant
plotResiduals(dharm_at, train$IABEF_2000*train$NDVIBEF_2000)# .25 is significant
#We have a problem, for many covariates one quantile is significant

hist(dharm_at)
#looks pretty flat >> thumbs up!

##Test for temporal autocorrelation
dharm_at_auto = recalculateResiduals(dharm_at, group = train$Fecha)
testTemporalAutocorrelation(dharm_at_auto, time =  unique(train$Fecha))
#doesn't  seem to be a problem

#Spatial Autocorrelation
dharm_at_spatial <- recalculateResiduals(dharm_at, group = train$Norte)
testSpatialAutocorrelation(dharm_at_spatial, 
                           x =  aggregate(train$Oeste, list(train$trap), mean)$x, 
                           y = aggregate(train$Norte, list(train$trap), mean)$x)
#No, spatial autocorrelation