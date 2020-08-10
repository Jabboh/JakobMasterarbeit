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
#b. Correlations between the responses / Residual Correlation Parameter of gjam
#c. Response Curves
#d. Variable importance
#6.2. Out-of-Sample
#a. Conditional Predictions of gjam vs. Unconditional Predictions of 
#gjam vs. predictions of univariate model
#b. Comparing the uncertainty of the different "prediction"-types
rm(list=ls())
setwd("C:\\Users\\jakob\\Documents\\JakobMasterarbeit\\Data")
#install.packages("clusteval")

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
                     I(NDVIBEF_2000^2), data = train, refresh = 0,
                   family = binomial(link = "probit"),init_r = .7, seed = 333)

#b. fitting the model for Anopheles
fit_at <- stan_glm(Anatropre ~ (IA_500 + NDVI_500 + IABEF_2000 + NDVIBEF_2000)^2 + 
                     Mes + I(IA_500^2) + I(NDVI_500^2) + I(IABEF_2000^2) +
                     I(NDVIBEF_2000^2), data = train, refresh = 0,
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
plotResiduals(dharm_at, train$IA_500) #.25 quantile is significant + combined quantile
#test also significant
plotResiduals(dharm_at, train$NDVI_500)#no significant problems 
plotResiduals(dharm_at, train$IABEF_2000)#significant problem(.25 quantile), but combined 
#quantiles insignificant
plotResiduals(dharm_at, train$NDVIBEF_2000)#.25 quantile is significant + combined quantile
#test also significant
plotResiduals(dharm_at, train$Mes) #looks good
#the quadratic forms
plotResiduals(dharm_at, (train$IA_500)^2) # .25 quantile is significant, but combined 
#quantile test insignificant
plotResiduals(dharm_at, (train$NDVI_500)^2) # .25 quantile is significant + comined
#significant
plotResiduals(dharm_at, (train$IABEF_2000)^2) # .25 quantile is significant, but combined 
#quantile test insignificant
plotResiduals(dharm_at, (train$NDVIBEF_2000)^2) # .25 quantile is significant + combined
#significant
#the interactions
plotResiduals(dharm_at, train$IA_500*train$NDVI_500)#.25 quantile is significant + 
#combined significant
plotResiduals(dharm_at, train$IA_500*train$IABEF_2000) # .25 quantile is significant
# + combined significant
plotResiduals(dharm_at, train$IA_500*train$NDVIBEF_2000) # .25 quantile is significant
#+ combined significant
plotResiduals(dharm_at, train$NDVI_500*train$IABEF_2000) # .25 and .5 quantile is significant
#but combined quantile test insignificant
plotResiduals(dharm_at, train$NDVI_500*train$NDVIBEF_2000) # .25 is significant + 
#combined significant
plotResiduals(dharm_at, train$IABEF_2000*train$NDVIBEF_2000)# .25 is significant +
#combined significant
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

#c.gjam DHARMa-analysis

#Prepare a newdata object for the gjam-posterior predictions



#I need to feed the Interaction and quadratic terms directly to the predict function
#So, I have to calculate these covariates first and store them as separate columns
train_gjdh <- as.data.frame(poly(as.matrix(train[,c("IA_500", "NDVI_500", "IABEF_2000",
                                                    "NDVIBEF_2000")]), degree = 2, raw = T))

#giving it the appropriate names
names(train_gjdh) <- c("IA500", "I(IA500^2)", "NDVI500", "IA500:NDVI500",
                       "I(NDVI500^2)", "IABEF2000", "IA500:IABEF2000",
                       "NDVI500:IABEF2000", "I(IABEF2000^2)", "NDVIBEF2000",
                       "IA500:NDVIBEF2000", "NDVI500:NDVIBEF2000", 
                       "IABEF2000:NDVIBEF2000", "I(NDVIBEF2000^2)")
#combining it with the other covariates
train_gjdh <- cbind(train_gj, train_gjdh)

#add intercept to train data (for some reason we need to this, otherwise gjamPredict
#won't accept it as xdata)
train_gjdh$intercept <- rep(1 , nrow(train)) 

# getting rid of all the unused variables
train_gjdh <- train_gjdh[,c("intercept", "IA500", "I(IA500^2)", "NDVI500", "IA500:NDVI500",
                            "I(NDVI500^2)", "IABEF2000", "IA500:IABEF2000",
                            "NDVI500:IABEF2000", "I(IABEF2000^2)", "NDVIBEF2000",
                            "IA500:NDVIBEF2000", "NDVI500:NDVIBEF2000", 
                            "IABEF2000:NDVIBEF2000", "I(NDVIBEF2000^2)", "Mayo", "Junio", "Julio",
                            "Agosto", "Septiembre", "Fecha", "trap", "Norte", "Oeste")] 

#I need to feed the Interaction and quadratic terms directly to the predict function
#So, I have to calculate these covariates first and store them as separate columns

#define the modelling setting
newdata <- list(xdata = train_gjdh, nsim = 4000)
#calculating the in-sample predictions (simulations)
sim <- gjamPredict(output = joint, newdata = newdata)

# The single draws from the predictive posterior are stored in ychains.
# It gives us values in terms of probabilities (the CDF right of the origin).
#In order to obtain simulated y-values on the observational scale that include the 
#residual error, we must draw from a binomial distribution with the probility that we
#drew from the predictive posterior distribution (rbinom does the job for us)

sims <- matrix(rbinom(n = length(sim$ychains),size = 1, sim$ychains),
               nrow = nrow(sim$ychains)) 

#creating the DHARMa object
dharm_gj_un <- createDHARMa(simulatedResponse = t(sims),
                            observedResponse = append(fit_cp$y, fit_at$y, after = length(fit_cp$y)),
                            fittedPredictedResponse = append(sim$sdList$yMu[,1], sim$sdList$yMu[,2],
                                                             after = length(sim$sdList$yMu[,1])),
                            integerResponse = T, seed = 333,method = "PIT")

plot(dharm_gj_un)

#1. QQ-plot looks good >> scaled residuals follow the expected uniform distribution
#2. KS-test also supports our assumption that the scaled residuals are uniformly 
#distributed >> our model is correctly specified.
#3. Residual vs. predicted quantile deviations (significant): The scaled residuals seem to increase
#with the predictions.

#Plotting the residuals against all covariates to check whether we specified the
# functional relationships correctly.

#Plot residuals against all covariates
plotResiduals(dharm_gj_un, rep(train$IA_500,2)) #looks ok, no significant problems
plotResiduals(dharm_gj_un, rep(train$NDVI_500, 2))#significant problem with the .5 quantile
#but combined quantile test is not signfificant
plotResiduals(dharm_gj_un, rep(train$IABEF_2000, 2))#no significant problems
plotResiduals(dharm_gj_un, rep(train$NDVIBEF_2000, 2))#significant problems for all
#quantiles and the combined test
plotResiduals(dharm_gj_un, rep(train$Mes, 2)) #looks good
#the quadratic forms
plotResiduals(dharm_gj_un, rep((train$IA_500)^2, 2)) # looks ok
plotResiduals(dharm_gj_un, rep((train$NDVI_500)^2, 2)) # looks ok
plotResiduals(dharm_gj_un, rep((train$IABEF_2000)^2, 2)) # looks ok
plotResiduals(dharm_gj_un, rep((train$NDVIBEF_2000)^2, 2)) # looks ok
#the interactions
plotResiduals(dharm_gj_un, rep(train$IA_500*train$NDVI_500, 2))#looks ok
plotResiduals(dharm_gj_un, rep(train$IA_500*train$IABEF_2000, 2)) # looks ok
plotResiduals(dharm_gj_un, rep(train$IA_500*train$NDVIBEF_2000, 2)) # looks ok
plotResiduals(dharm_gj_un, rep(train$NDVI_500*train$IABEF_2000, 2)) # looks ok
plotResiduals(dharm_gj_un, rep(train$NDVI_500*train$NDVIBEF_2000, 2)) # looks ok
plotResiduals(dharm_gj_un, rep(train$IABEF_2000*train$NDVIBEF_2000, 2))# no problems
#Overall, only two plots had significant problems and for only one of them, the
#combined test was significant . Due to the high number of plots,
#this is somewhat  expected (we performed 15*3 = 45 significance tests >> 
#hence we would expect roughly 2 significant plots (95%-level)) >> So, we should be
#good! :)

hist(dharm_gj_un)
#looks pretty flat >> thumbs up!

##Test for temporal autocorrelation
dharm_gj_un_auto = recalculateResiduals(dharm_gj_un, group = rep(train_gjdh$Fecha, 2))
testTemporalAutocorrelation(dharm_gj_un_auto, time =  unique(train_gjdh$Fecha))
#doesn't  seem to be a problem

#Spatial Autocorrelation

#I do not have the coordinates for trap "M29" --> one NA in the coordinates >> We need to remove that row
dharm_gj_un_spatial <- recalculateResiduals(dharm_gj_un,
                                            group = rep(train_gjdh$Norte, 2))
testSpatialAutocorrelation(dharm_gj_un_spatial, 
                           x =  aggregate(train$Oeste, list(train$trap), mean)$x, 
                           y = aggregate(train$Norte, list(train$trap), mean)$x)
#No spatial autocorrelation >> Yeah!

#4. Finding the "best" model
#I proceed as follows. I drop each environmental covariate(month, NDVI, IA...) 
#and the related terms (quadratic + interactions) separately and choose the model with
#the lowest WAIC (the best reduced model). If this WAIC is higher (meaning worse)
#than the WAIC of the complex model, I repeat the steps until the more complex model has
#a lower WAIC (meaning better model).
#At the end of this first step,
#I have chosen the principal covariates of my model. Next, I want to choose the
#functional relationship of them with regard to the response. For this, I first
#check whether dropping the interaction terms reduces the WAIC. So I drop them
#iteratevly and again choose the model with the lowest AIC. Afterwards I do 
#the same for the quadratic forms of the covariates.
#x: fitted model for which coefficients should be dropped
#y: names of the covariates 
#Loading the step functions I created
source("WaicStep.R")

#Finding the "best" model for Culex perexiguus: Dropping environmental 
#covariates and all its associations

#defining a vector with all the names of the covariates
z <- c("IA_500", "NDVI_500", "IABEF_2000", "NDVIBEF_2000", "Mes")

step_waic_entCov(fit_cp, z)
#The function dropped two variables NDVI_500 and IABEF_2000. So, the new model
#includes IA_500, NDVIBEF_2000 and the interactions as well as the quadratic terms
#and the factor mes.
#Running this model again and storing it (we need it as an input for the
#interquad function):
fitr_cp <-  stan_glm(Cxperpre ~ (IA_500 + NDVIBEF_2000)^2 + Mes +
                    I(IA_500^2) + I(NDVIBEF_2000^2), data = train,
                    refresh = 0, family = binomial(link = "probit"), 
                    init_r = .7, seed = 333)

#defining a names vector
z <- c("IA_500", "NDVIBEF_2000")
interquad(fitr_cp, z)
#So, we drop the only remaining interaction term, but keep all the quadratic
#terms. Hence the new CP model, after reducing its complexity on the basis
#of WAIC:
fitp_cp <- stan_glm(Cxperpre ~ Mes + IA_500 + NDVIBEF_2000 +
                    I(IA_500^2) + I(NDVIBEF_2000^2), data = train,
                    family = binomial(link = "probit"), init_r = .7, 
                    seed = 333)
  
#Finding the "best" model for Anopheles troparvus: Dropping environmental 
#covariates and all its associations

#defining a vector with all the names of the covariates
z <- c("IA_500", "NDVI_500", "IABEF_2000", "NDVIBEF_2000", "Mes")

step_waic_entCov(fit_at, z)
#The function dropped two variables NDVI_500 and IABEF_2000. So, the new model
#includes IA_500, NDVIBEF_2000 and the interactions as well as the quadratic terms
#and the factor mes. (Same model as for CP)

#Running this model again and storing it (we need it as an input for the
#interquad function):
fitr_at <-  stan_glm(Anatropre ~ (IA_500 + NDVIBEF_2000)^2 + Mes +
                       I(IA_500^2) + I(NDVIBEF_2000^2), data = train,
                     refresh = 0, family = binomial(link = "probit"), 
                     init_r = .7, seed = 333)

#defining a names vector
z <- c("IA_500", "NDVIBEF_2000")
interquad(fitr_at, z)
#So, we drop the only remaining interaction term and the quadratic
#term of IA_500. 
fitp_at <- stan_glm(Anatropre ~ IA_500 + NDVIBEF_2000 + Mes +
                    I(NDVIBEF_2000^2), data = train,
                    refresh = 0, family = binomial(link = "probit"), 
                    init_r = .7, seed = 333)
#the "best" model specification according to CP, yields the following model
#for at.

fitcp_at <- stan_glm(Anatropre ~ Mes + IA_500 + NDVIBEF_2000 +
                       I(IA_500^2) +  I(NDVIBEF_2000^2), data = train,
                    refresh = 0, family = binomial(link = "probit"), 
                    init_r = .7, seed = 333)
#comparing their WAICs
waic(fitcp_at)
waic(fitp_at)
#the difference in WAIC is 1

#the "best" model specification according to AT, yields the follwoing model
#for CP
fitat_cp <- stan_glm(Cxperpre ~ IA_500 + NDVIBEF_2000 + Mes +
                      I(NDVIBEF_2000^2), data = train,
                    refresh = 0, family = binomial(link = "probit"), 
                    init_r = .7, seed = 333)
waic(fitat_cp)
waic(fitp_cp)
#The difference in WAIC is 1.3

#So, the "best" models for at and cp only have one difference: the quadratic 
#term of IA_500. Since we need the same models to compare the results, we
#need to choose one of the two models for both species. There are two reasons
#for choosing the model with the quadratic term: 1. Including more rather
#than fewer variables, is a safer strategy with regard to misspecifications
#of the model and 2. the difference in WAIC between the two models is smaller
#for AT, meaning that including this variable for AT has a lesser negative 
#impact on WAIC than excluding the quadratic term has on the WAIC in the case
#of CP.

#Hence, our three final models are the following:

fit_fin_cp <- fitp_cp
fit_fin_at <- fitcp_at
#gjam
#Define the model settings

#define model/algorithm parameters: 4000 gibbs steps (equivalent to the MCMCs draws
#in the rstanarm model fitting) + burnin of 1000
ml   <- list(ng = 4000, burnin = 1000, typeNames = types)

#runnig GJAM
joint_fin <- gjam(~  Mayo + Junio + Julio + Agosto + Septiembre + IA_500 +
                NDVIBEF_2000 + I(IA_500^2) + I(NDVIBEF_2000^2),
                ydata = y_train, xdata = train_gj, modelList = ml)
joint_fin$fit$DIC
#The DIC (related to WAIC) is also lower for new gjam model. 
####5. Doing DHARMa on the final models
#For cp
dharm_fin_cp <- createDHARMa(simulatedResponse = t(posterior_predict(fit_fin_cp)), observedResponse = fit_fin_cp$y,
                         fittedPredictedResponse = posterior_predict(fit_fin_cp) %>% apply(2, median), integerResponse = T, seed = 333,
                         method = "PIT")
plot(dharm_fin_cp)
#1. QQ-plot looks good >> scaled residuals follow the expected uniform distribution
#2. KS-test also supports our assumption that the scaled residuals are uniformly 
#distributed >> our model is correctly specified.
#3. Boxplot: For both model predictions (absent or present), the scaled residuals are
#approximately equally distributed >> no bias

#Plotting the residuals against all covariates to check whether we specified the
# functional relationships correctly.

#Plot residuals against all covariates
plotResiduals(dharm_fin_cp, train$IA_500) #looks ok, no significant problems
plotResiduals(dharm_fin_cp, train$NDVIBEF_2000)#no problems
plotResiduals(dharm_fin_cp, train$Mes) #looks good
#the quadratic forms
plotResiduals(dharm_fin_cp, (train$IA_500)^2) # looks ok
plotResiduals(dharm_fin_cp, (train$NDVIBEF_2000)^2) # looks ok
#Overall, no plots had significant problems >> Can ride with our specifications

hist(dharm_fin_cp)
#looks pretty flat >> thumbs up!

##Test for temporal autocorrelation
dharm_fin_cp_auto = recalculateResiduals(dharm_fin_cp, group = train$Fecha)
testTemporalAutocorrelation(dharm_fin_cp_auto, time =  unique(train$Fecha))
#doesn't  seem to be a problem

#Spatial Autocorrelation

dharm_fin_cp_spatial <- recalculateResiduals(dharm_fin_cp, group = train$Norte)
testSpatialAutocorrelation(dharm_fin_cp_spatial, 
                           x =  aggregate(train$Oeste, list(train$trap), mean)$x, 
                           y = aggregate(train$Norte, list(train$trap), mean)$x)
#No spatial autocorrelation >> Yeah!
#Overall, DHARMa does not flag any issues with the model

#For at
dharm_fin_at <- createDHARMa(simulatedResponse = t(posterior_predict(fit_fin_at)), observedResponse = fit_fin_at$y,
                             fittedPredictedResponse = posterior_predict(fit_fin_at) %>% apply(2, median), integerResponse = T, seed = 333,
                             method = "PIT")
plot(dharm_fin_at)
#1. QQ-plot looks good >> scaled residuals follow the expected uniform distribution
#2. KS-test also supports our assumption that the scaled residuals are uniformly 
#distributed >> our model is correctly specified.
#3. Boxplot: For both model predictions (absent or present), the scaled residuals are
#approximately equally distributed >> no bias

#Plotting the residuals against all covariates to check whether we specified the
# functional relationships correctly.

#Plot residuals against all covariates
plotResiduals(dharm_fin_at, train$IA_500) #.25 quantile is significant, but
#combined quantile test insignificant
plotResiduals(dharm_fin_at, train$NDVIBEF_2000)#.25 quantile is significant
#and combined quantile test also significant
plotResiduals(dharm_fin_at, train$Mes) #looks good
#the quadratic forms
plotResiduals(dharm_fin_at, (train$IA_500)^2) #.25 quantile is significant
#and combined quantile test also significant
plotResiduals(dharm_fin_at, (train$NDVIBEF_2000)^2) #.25 quantile is significant
#and combined quantile test also significant

#Overall, significant problems with the predictors

hist(dharm_fin_at)
#looks pretty flat >> thumbs up!

##Test for temporal autocorrelation
dharm_fin_at_auto = recalculateResiduals(dharm_fin_at, group = train$Fecha)
testTemporalAutocorrelation(dharm_fin_at_auto, time =  unique(train$Fecha))
#doesn't  seem to be a problem

#Spatial Autocorrelation

dharm_fin_at_spatial <- recalculateResiduals(dharm_fin_at, group = train$Norte)
testSpatialAutocorrelation(dharm_fin_at_spatial, 
                           x =  aggregate(train$Oeste, list(train$trap), mean)$x, 
                           y = aggregate(train$Norte, list(train$trap), mean)$x)
#No spatial autocorrelation >> Yeah!
#Overall, DHARMa does not flag any issues with the model

#for gjam
#Prepare a newdata object for the gjam-posterior predictions

# getting rid of all the unused variables
train_gjdh <- train_gjdh[,c("intercept", "IA500", "I(IA500^2)", "NDVIBEF2000",
                             "I(NDVIBEF2000^2)", "Mayo", "Junio", "Julio",
                            "Agosto", "Septiembre", "Fecha", "trap", "Norte", "Oeste")] 

#define the modelling setting
newdata <- list(xdata = train_gjdh, nsim = 4000)
#calculating the in-sample predictions (simulations)
sim_fin <- gjamPredict(output = joint_fin, newdata = newdata)

# The single draws from the predictive posterior are stored in ychains.
# It gives us values in terms of probabilities (the CDF right of the origin).
#In order to obtain simulated y-values on the observational scale that include the 
#residual error, we must draw from a binomial distribution with the probility that we
#drew from the predictive posterior distribution (rbinom does the job for us)

sims_fin <- matrix(rbinom(n = length(sim_fin$ychains),size = 1, sim_fin$ychains),
               nrow = nrow(sim_fin$ychains)) 

#creating the DHARMa object
dharm_gj_un_fin <- createDHARMa(simulatedResponse = t(sims_fin),
                            observedResponse = append(fit_cp$y, fit_at$y, after = length(fit_cp$y)),
                            fittedPredictedResponse = append(sim_fin$sdList$yMu[,1], sim_fin$sdList$yMu[,2],
                                                             after = length(sim_fin$sdList$yMu[,1])),
                            integerResponse = T, seed = 333,method = "PIT")

plot(dharm_gj_un_fin)

#1. QQ-plot looks good >> scaled residuals follow the expected uniform distribution
#2. KS-test also supports our assumption that the scaled residuals are uniformly 
#distributed >> our model is correctly specified.
#3. Residual vs. predicted quantile deviations (significant): The scaled residuals seem to increase
#with the predictions.

#Plotting the residuals against all covariates to check whether we specified the
# functional relationships correctly.

#Plot residuals against all covariates
plotResiduals(dharm_gj_un_fin, rep(train$IA_500,2)) #looks ok, no significant problems
plotResiduals(dharm_gj_un_fin, rep(train$NDVIBEF_2000, 2))#significant problems for all
#quantiles and the combined test
plotResiduals(dharm_gj_un_fin, rep(train$Mes, 2)) #looks good
#the quadratic forms
plotResiduals(dharm_gj_un_fin, rep((train$IA_500)^2, 2)) # looks ok
plotResiduals(dharm_gj_un_fin, rep((train$NDVIBEF_2000)^2, 2)) # looks ok

#Overall, only one plot has significant problems. Some significant results
#are expected:)

hist(dharm_gj_un_fin)
#looks pretty flat >> thumbs up!

##Test for temporal autocorrelation
dharm_gj_un_fin_auto = recalculateResiduals(dharm_gj_un_fin, group = rep(train_gjdh$Fecha, 2))
testTemporalAutocorrelation(dharm_gj_un_fin_auto, time =  unique(train_gjdh$Fecha))
#doesn't  seem to be a problem

#Spatial Autocorrelation
dharm_gj_un_fin_spatial <- recalculateResiduals(dharm_gj_un_fin,
                                            group = rep(train_gjdh$Norte, 2))
testSpatialAutocorrelation(dharm_gj_un_fin_spatial, 
                           x =  aggregate(train$Oeste, list(train$trap), mean)$x, 
                           y = aggregate(train$Norte, list(train$trap), mean)$x)
#No spatial autocorrelation >> Yeah!

####6. Results
###6.1. In-sample:
##a. Comparing the coefficients: Size, Significance and credibility intervals
#Plotting all the coefficients in one plot per species
#CP
#https://github.com/stan-dev/bayesplot/issues/232
#The posteriors of both models
posterior_1 <- as.matrix(fit_fin_cp) #rstanarm nimmt den median als Schätzer für den Parameter
posterior_2 <- as.matrix(joint_fin$chains$bgibbsUn[,1:10])
colnames(posterior_2) <- colnames(posterior_1)
combined <- rbind(mcmc_intervals_data(posterior_1, prob_outer = .95), mcmc_intervals_data(posterior_2, prob_outer = .95))
combined$model <- rep(c("Univariate Model", "Multivariate Model"), each = ncol(posterior_1))
#prob_outer defines the credibility interval in our plot (here .975)
theme_set(bayesplot::theme_default())
pos <- position_nudge(y = ifelse(combined$model == "Multivariate Model", 0, 0.1))
coef_cp <- ggplot(combined, aes(x = m, y = parameter, color = model)) + 
  geom_point(position = pos) +
  geom_vline(xintercept = 0, linetype="dotted", color = "black", size=.5) +
  geom_errorbar(aes(xmin = ll, xmax = hh), position = pos, width = .1) +
  ggtitle("Univariate vs. Multivariate Coefficients and their 95 % - Credibility Intervals \n for Culex Perexiguus") + 
  xlab("Value") + ylab("Coefficient") + labs(color="Model") 
#They look pretty much the same 

#for at
posterior_1 <- as.matrix(fit_fin_at) #rstanarm nimmt den median als Schätzer für den Parameter
posterior_2 <- as.matrix(joint_fin$chains$bgibbsUn[,-(1:10)])
colnames(posterior_2) <- colnames(posterior_1)
combined <- rbind(mcmc_intervals_data(posterior_1, prob_outer = .95), mcmc_intervals_data(posterior_2, prob_outer = .95))
combined$model <- rep(c("Univariate Model", "Multivariate Model"), each = ncol(posterior_1))
#prob_outer defines the credibility interval in our plot (here .975)
pos <- position_nudge(y = ifelse(combined$model == "Multivariate Model", 0, 0.1))
coef_at <- ggplot(combined, aes(x = m, y = parameter, color = model)) + 
  geom_point(position = pos) +
  geom_errorbar(aes(xmin = ll, xmax = hh), position = pos, width = .1) +
  geom_vline(xintercept = 0, linetype="dotted", color = "black", size=.5) +
  ggtitle("Univariate vs. Multivariate Coefficients and their 95 % - Credibility Intervals \n for Culex Perexiguus") + 
  xlab("Value") + ylab("Coefficient") + labs(color="Model") 
#They look pretty much the same 

##the hard numbers: Tables with the coefficients and SEs
#Significance: (guck eigentlich am besten deinen Plot an!)
summary(joint_fin)
summary(fit_fin_cp)
summary(fit_fin_at)

#Comparison of the coefficients of CP
#Make a table (dataframe) to store different coefficients and SEs according to the models
q = length(fit_fin_cp$coefficients) # Number of predictors
#For Culex Perexiguus
#For the coefficients
cof_sum_cp <- data.frame(matrix(ncol = 2, nrow = q))
colnames(cof_sum_cp) <- c("Coefficients_sin", "Coefficients_gjam")
rownames(cof_sum_cp) <- names(fit_fin_cp$coefficients)

#For the SEs
se_sum_cp <- data.frame(matrix(ncol = 2, nrow = q))
colnames(se_sum_cp) <- c("SE_sin", "SE_gjam")
rownames(se_sum_cp) <- names(fit_fin_cp$coefficients)

###For Anopheles troparvus
#For the coefficients
cof_sum_at <- cof_sum_cp
#for the SEs
se_sum_at <- se_sum_cp

#Filling the tables accordingly
cof_sum_cp$Coefficients_sin <-  fit_fin_cp$coefficients
cof_sum_cp$Coefficients_gjam <- joint_fin$parameters$betaMu[,"Cxperpre"]
cof_sum_at$Coefficients_sin <-  fit_fin_at$coefficients
cof_sum_at$Coefficients_gjam <- joint_fin$parameters$betaMu[,"Anatropre"]

se_sum_cp$SE_sin <-  fit_fin_cp$ses
se_sum_cp$SE_gjam <- joint_fin$parameters$betaSe[,"Cxperpre"]
se_sum_at$SE_sin <-  fit_fin_at$ses
se_sum_at$SE_gjam <- joint_fin$parameters$betaSe[,"Anatropre"]

#Coefficients and SEs for Culex perexiguus
cof_sum_cp
se_sum_cp
#Everything looks pretty similar, as we expected (We expexted the environmental coefficients to
#be the same. Difference between coefficients way smaller than according SEs. 
#You can also check that by looking at the graph

#Coefficients and SEs for Anopheles troparvus
cof_sum_at
se_sum_at
#Everything looks pretty similar as expected. Difference between coefficients way smaller than according SEs.
#only exception is september dummy, but for that the SE is also very large.

###Correlation between the responses (raw vs. residual correlation)
#raw correlation
cor(spec$Cxperpre, spec$Anatropre) #pearson correlation coefficient

#a better measure for similarity of binary data (Jaccard Index 
#(Intersection(A,B)/Union(A,B)) )
library('clusteval')
cluster_similarity(spec$Cxperpre,spec$Anatropre, similarity="jaccard",
                   method="independence")
#not a big difference compared to raw correlation

#Residual Correlations
joint_fin$parameters$corMu #It's the pearson correlation coefficient (for two
#continuous variables)

#residual correlation is smaller >> reason: shared environmental response is
#accounted for (Most coefficient have the same sign (only a couple of month
#dummies are different)! But there is still a lot of unexplained Correlation
#(are you sure? What's the natural level? I guess 0)

###c. Response Curves
