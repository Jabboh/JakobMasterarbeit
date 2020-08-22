################################Analysis of Culex perexiguus & Anopheles troparvus:
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

#Selecting our two species (perexiguus & Anopheles troparvus) for the analysis. Our
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

####Fitting the most complex models
#The covariates Inundation Area and NDVI are available on different spatial and
#temporal scales. For the temporal scales, I include both scales (current +
#month before) in the analyis. For the spatial scale, I select one of the five
#alternatives (100 m buffer, 250 m, 500 m, 1000 m and 2000m). I decided which scale
#to choose on the basis of the models in Roiz's paper. He already has "the most
#appropriate" specifications for his Generalized Linear Mixed Models (GLMM) with a 
#binomial error distribution and binomial link function for Culex perexiguus & 
#Anopheles troparvus. So, for example since he used "NDVI_500" for his Culex perexiguus
#model, I use also the "NDVI_500" in all my models. Justifications for my choice of 
#spatial scales:
#IA_500: In Roiz's Culex perexiguus model
#NDVI_500: In Roiz's Culex perexiguus model
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

#a. fitting the model for Culex perexiguus
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
#a. Culex perexiguus
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
                    family = binomial(link = "probit"), init_r = 1.4, 
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
                    init_r = 1.4, seed = 333)
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
  ggtitle("Univariate vs. Multivariate Coefficients and their 95 % - Credibility Intervals \n for Culex perexiguus") + 
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
  ggtitle("Univariate vs. Multivariate Coefficients and their 95 % - Credibility Intervals \n for Culex perexiguus") + 
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
#For Culex perexiguus
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
##For Culex perexiguus

###Response Curve für IA_500 :predict the probability of presence for the different ia-values,
#holding the other covariates at their mean (0) or at a constant factor level (mes)

#a sequence of IA_500-values from min to max with 50 steps
ia <- seq(min(df$IA_500), max(df$IA_500), length.out = 50)
#names of covariates
nc <- c("IA_500", "NDVIBEF_2000", "Mes")
#the mean of the covariates is 0, bc I normalized them
#creating  x-data dataframe with all the avareges (which are zero)
data <- as.data.frame(matrix(0, ncol = length(nc), nrow = length(ia)))
names(data) <- nc  
#replace IA_500 0s with the sequence from above
data$IA_500 <-ia
#Next, we want to create a dataframe with "data" for all the different months (bc we want to have
#a response curve for every month separately)
#initialize empty data frame that will be filled in the loop over the months
xdata <- as.data.frame(matrix(, ncol = length(nc), nrow = 0))
names(xdata) <- nc  
#run a loop over all the months and append month after month (with the same 
#covariates) to the dataframe xdata
for (i in levels(df$Mes)){
  d<- data
  d$Mes <- i
  xdata <- rbind(xdata, d)
}
#Converting Mes to a factor variable
xdata$Mes <- factor(xdata$Mes, levels =c("Abril", "Mayo", "Junio", "Julio", "Agosto",
                                         "Septiembre"))
#Run the posterior simulations with xdata, we want the predictions on the "probability scale"
#That's why we use posterior_linpred with transform = T
d_res <- posterior_linpred(fit_fin_cp, transform = T, newdata = xdata, seed = 333)
#getting the predictions on the "probability scale, by taking the mean per 
#observation/column
univariate <- colMeans(d_res)
#adding these predictions to data frame
ggd <- cbind(xdata, univariate)
#maiking the plot
response_cp_uv <- ggplot(data = ggd, aes(x = IA_500, y = univariate, color = Mes)) +
  geom_point() + ylab("Predicted Probability") + xlab("IA_500 in Standard Units") +
  ggtitle ("Response Curve of Inundation Area for Culex perexiguus in Univariate Model") +
  guides(color=guide_legend(title="Month"))

##Doing the same thing for gjam

#for unconditional predictions
#Converting the Mes-Variable to a bunch of dummies as demanded by gjam
source("gjamDummy.R")

dr_gj <- gjamDummy(xdata$Mes, xdata)
#rename certain columns to coefficient names of joint_fin output
names(dr_gj)[1:2] <- c("IA500", "NDVIBEF2000")
#add an intercept to the data frame
dr_gj$intercept <- rep(1 , nrow(xdata)) 
#add already calculated quadratic terms
dr_gj$"I(IA500^2)" <- (dr_gj$IA500)^2
dr_gj$"I(NDVIBEF2000^2)" <- 0
#add a superfluous line of 1s (gjamPredict doesnt run, if there is no
#variation in the covariate) >> we will delete this row after making the 
#predictions
dr_gj <- rbind(dr_gj,1)
#delete Mes variable (also needed for gjamPredict to run)
dr_gj$Mes <- NULL
#define the modelling settings for posterior simulations with dr_gj
newdata <- list(xdata = dr_gj, nsim = 4000)
#calculating the in-sample predictions (simulations)
sim <- gjamPredict(output = joint_fin, newdata = newdata)

#take the first 300 rows as predictions (deleting the last row); Predictions are the means of the
#ychains
pred_gj_un <- sim$sdList$yMu[1:300,1] 
dr_gj <- dr_gj[-nrow(dr_gj),]
#readding the Mes column
dr_gj$Mes <- ggd$Mes
#making the ggplot dataframe
ggd_gj_un <- cbind(dr_gj, pred_gj_un)
response_cp_un <- ggplot(data = ggd_gj_un, aes(x = IA500, y = pred_gj_un, color = Mes)) +
  geom_point() + ylab("Predicted Probability") + xlab("IA_500 in Standard Units") +
  ggtitle ("Response Curve of Inundation Area for Culex perexiguus in GJAM") +
  guides(color=guide_legend(title="Month"))

##Plot gjam and rstanarm response curves in the same plot

#prepare the dataframe to feed into ggplot (just add predictions of gjam to dataframe of rstanarm)

ggd$multivariate <- ggd_gj_un$pred_gj_un

#Do the ggplot 

response_cp_uv_un <- ggplot(data = ggd, aes(x = IA_500, color = Mes)) +
  geom_point(aes(y = univariate, shape = "univariate")) + 
  geom_point(aes(y = multivariate, shape = "multivariate")) + ylab("Predicted Probability") + 
  xlab("IA_500 in Standard Units") +
  ggtitle ("Univariate vs. Unconditional Multivariate Response Curves of Inundation Area for Culex perexiguus") +
  guides(color=guide_legend(title="Month"), shape = guide_legend(title="Model"))

#Why are the response curves different? Does that contradict our hypothesis that unconditional
#predictions dont differ from univariate predictions? Why do they cross? (heterogeneous effect (
#of month depends on value of IA_500)) Why is for september multivariate probability more positive,
#but for all the other months the univariate model

#DOing the same thing with conditional predictions of gjam

#replicate dr_gj: We want to store both conditional prediction in one data frame (so we run gjamPredict
#only once and obtain both conditional predictions
dr_gj_con <- rbind(dr_gj, dr_gj)


#delete Mes variable (also needed for gjamPredict to run)
dr_gj_con$Mes <- NULL
#add a superfluous line of 1s (gjamPredict doesnt run, if there is no
#variation in the covariate) >> we will delete this row after making the 
#predictions
dr_gj_con <- rbind(dr_gj_con,1)
#create the PA-values for the Anopheles, the species we condition on (we create an addional one 
#so that it has the same length as dr_gj_con)
ydata <- tibble(c(rep(0, 300), rep(1, 301)))
names(ydata) <- "Anatropre"


#define the modelling settings for posterior simulations with the conditional species as present
newdata <- list(xdata = dr_gj_con, ydataCond = ydata, nsim = 4000)

#calculating the in-sample predictions (simulations)
sim <- gjamPredict(output = joint_fin, newdata = newdata)

#take the first 600 rows as predictions (deleting the last row); Predictions are the means of the
#ychains
pred_gj_con <- sim$sdList$yMu[1:600,1] 

#adding the conditional presence of Anotrophes to our data storage data frame
dr_gj_con$cond <- factor(ydata$Anatropre, labels = c("absent", "present"))
#deleting the last row of our data (remember: we added only bc this way the predict function runs)
dr_gj_con <- dr_gj_con[-601,]
#readding the Mes column
dr_gj_con$Mes <- ggd$Mes

dr_gj_con <- cbind(dr_gj_con, pred_gj_con)
response_cp_gj_con <- ggplot(data = dr_gj_con, aes(x = IA500, y = pred_gj_con, color = Mes, shape = cond)) +
  geom_point() + ylab("Predicted Probability") + xlab("IA_500 in Standard Units") +
  ggtitle ("Response Curve of Inundation Area for Culex perexiguus in GJAM \n Conditional on Presence of Anophles Atroparvus") +
  guides(color=guide_legend(title="Month"), shape = guide_legend(title = "Anopheles atroparvus"))

#strange results, ordering of the effect of month on pred, changes depending on whether AT is present
#or not. Does that make sense? Is that even possible? >> I actually do not think so (BUt only september
#is different!)

###Doing a plot of gjam conditional predictions and univariate prediction (We just take any month)

#taking June


#Preparing the input data
plots_uv_con <- vector(mode = "list", length = length(levels(xdata$Mes)))
for( i in levels(xdata$Mes)){
  #index variable
  j <-match(i, levels(xdata$Mes))
  nu <- nrow(xdata[xdata["Mes"] == i,])
  d <- rbind(xdata[xdata["Mes"] == i,], xdata[xdata["Mes"] == i,], xdata[xdata["Mes"] == i,])
  #creating one column that has the presence-absence of at for the gjam models and for the univariate
  #model the variable takes on another factor level (2)
  d$mode <- factor(c(rep(0, nu), rep(1, nu), rep(2, nu)), labels = c("GJAM with absent Anopheles", 
                                                                     "GJAM with present Anopheles",
                                                                     "Univariate Model"))
  #adding the predictions
  d$pred <- c(dr_gj_con$pred_gj_con[dr_gj_con$Mes == i], ggd$univariate[ggd$Mes == i])
  
  #Doing the ggplot
  plot <- ggplot(data = d, aes(x = IA_500, y = pred, shape = mode)) +
    geom_point() + ylab("Predicted Probability") + xlab("IA_500 in Standard Units") +
    ggtitle (paste0("Response Curve of Inundation Area for Culex perexiguus in GJAM \n Conditional on Presence of Anophles Atroparvus vs. Univariate Model for Month ", i)) +
    theme(legend.title = element_blank())
  plots_uv_con[[j]] <- plot
}
plots_uv_con[[2]]

###################################################
#for NDVIBEF#

#a sequence of NDVIBEF-values from min to max with 50 steps
nd <- seq(min(df$NDVIBEF_2000), max(df$NDVIBEF_2000), length.out = 50)
#names of covariates
nc <- c("IA_500", "NDVIBEF_2000", "Mes")
#the mean of the covariates is 0, bc I normalized them
#creating  x-data dataframe with all the avareges (which are zero)
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
for (i in levels(df$Mes)){
  d<- data
  d$Mes <- i
  xdata <- rbind(xdata, d)
}
#Converting Mes to a factor variable
xdata$Mes <- factor(xdata$Mes, levels =c("Abril", "Mayo", "Junio", "Julio", "Agosto",
                                         "Septiembre"))
#Run the posterior simulations with xdata, we want the predictions on the "probability scale"
#That's why we use posterior_linpred with transform = T
d_res <- posterior_linpred(fit_fin_cp, transform = T, newdata = xdata, seed = 333)
#getting the predictions on the "probability scale, by taking the mean per 
#observation/column
univariate <- colMeans(d_res)
#adding these predictions to data frame
ggd <- cbind(xdata, univariate)
#maiking the plot
response_cp_uv <- ggplot(data = ggd, aes(x = NDVIBEF_2000, y = univariate, color = Mes)) +
  geom_point() + ylab("Predicted Probability") + xlab("NDVIBEF_2000 in Standard Units") +
  ggtitle ("Response Curve of NDVI Before for Culex perexiguus in Univariate Model") +
  guides(color=guide_legend(title="Month"))

#Interesting overall shape >> kinda contradicts the niche concept >> we have maxima at the edges
##Doing the same thing for gjam

#for unconditional predictions
#Converting the Mes-Variable to a bunch of dummies as demanded by gjam and storing whole data in
#dr_gj
dr_gj <- gjamDummy(xdata$Mes, xdata)
#rename certain columns to coefficient names of joint_fin output
names(dr_gj)[1:2] <- c("IA500", "NDVIBEF2000")
#add an intercept to the data frame
dr_gj$intercept <- rep(1 , nrow(xdata)) 
#add quadratic terms
dr_gj$"I(IA500^2)" <- 0
dr_gj$"I(NDVIBEF2000^2)" <- (dr_gj$NDVIBEF2000)^2
#add a superfluous line of 1s (gjamPredict doesnt run, if there is no
#variation in the covariate) >> we will delete this row after making the 
#predictions
dr_gj <- rbind(dr_gj,1)
#delete Mes variable (also needed for gjamPredict to run)
dr_gj$Mes <- NULL
#define the modelling settings for posterior simulations with dr_gj
newdata <- list(xdata = dr_gj, nsim = 4000)
#calculating the in-sample predictions (simulations)
sim <- gjamPredict(output = joint_fin, newdata = newdata)

#take the first 300 rows as predictions (deleting the last row); Predictions are the means of the
#ychains
pred_gj_un <- sim$sdList$yMu[1:300,1] 
dr_gj <- dr_gj[-nrow(dr_gj),]
#readding the Mes column
dr_gj$Mes <- ggd$Mes
#making the ggplot dataframe
ggd_gj_un <- cbind(dr_gj, pred_gj_un)
response_cp_un <- ggplot(data = ggd_gj_un, aes(x = NDVIBEF2000, y = pred_gj_un, color = Mes)) +
  geom_point() + ylab("Predicted Probability") + xlab("NDVIBEF_2000 in Standard Units") +
  ggtitle ("Response Curve of NDVI Before for Culex perexiguus in GJAM") +
  guides(color=guide_legend(title="Month"))

##Plot gjam and rstanarm response curves in the same plot

#prepare the dataframe to feed into ggplot (just add predictions of gjam to dataframe of rstanarm)

ggd$multivariate <- ggd_gj_un$pred_gj_un

#Do the ggplot 

response_cp_uv_un <- ggplot(data = ggd, aes(x = NDVIBEF_2000, color = Mes)) +
  geom_point(aes(y = univariate, shape = "univariate")) + 
  geom_point(aes(y = multivariate, shape = "multivariate")) + ylab("Predicted Probability") + 
  xlab("NDVIBEF_2000 in Standard Units") +
  ggtitle ("Univariate vs. Multivariate Response Curves of NDVI Before for Culex perexiguus") +
  guides(color=guide_legend(title="Month"), shape = guide_legend(title="Model"))

#They have a whole different shape >> I dont't think they should behave differently, but apparently
#they do o_o



#DOing the same thing with conditional predictions of gjam

#replicate dr_gj: We want to store both conditional prediction in one data frame (so we run gjamPredict
#only once and obtain both conditional predictions
dr_gj_con <- rbind(dr_gj, dr_gj)


#delete Mes variable (also needed for gjamPredict to run)
dr_gj_con$Mes <- NULL
#add a superfluous line of 1s (gjamPredict doesnt run, if there is no
#variation in the covariate) >> we will delete this row after making the 
#predictions
dr_gj_con <- rbind(dr_gj_con,1)
#create the PA-values for the Anopheles, the species we condition on (we create an addional one 
#so that it has the same length as dr_gj_con)
ydata <- tibble(c(rep(0, 300), rep(1, 301)))
names(ydata) <- "Anatropre"


#define the modelling settings for posterior simulations with the conditional species as present
newdata <- list(xdata = dr_gj_con, ydataCond = ydata, nsim = 4000)

#calculating the in-sample predictions (simulations)
sim <- gjamPredict(output = joint_fin, newdata = newdata)

#take the first 600 rows as predictions (deleting the last row); Predictions are the means of the
#ychains
pred_gj_con <- sim$sdList$yMu[1:600,1] 

#adding the conditional presence of Anotrophes to our data storage data frame
dr_gj_con$cond <- factor(ydata$Anatropre, labels = c("absent", "present"))
#deleting the last row of our data (remember: we added only bc this way the predict function runs)
dr_gj_con <- dr_gj_con[-601,]
#readding the Mes column
dr_gj_con$Mes <- ggd$Mes

dr_gj_con <- cbind(dr_gj_con, pred_gj_con)
response_cp_gj_con <- ggplot(data = dr_gj_con, aes(x = NDVIBEF2000, y = pred_gj_con, color = Mes, shape = cond)) +
  geom_point() + ylab("Predicted Probability") + xlab("NDVIBEF_2000 in Standard Units") +
  ggtitle ("Response Curve of NDVI Before for Culex perexiguus in GJAM \n Conditional on Presence of Anophles Atroparvus") +
  guides(color=guide_legend(title="Month"), shape = guide_legend(title = "Anopheles atroparvus"))

#strange results, ordering of the effect of month on pred, changes depending on whether AT is present
#or not. Does that make sense? Is that even possible? >> I actually do not think so (BUt only september
#is different!)

###Doing a plot of gjam conditional predictions and univariate prediction (We just take any month)


#Preparing the input data
plots_uv_con <- vector(mode = "list", length = length(levels(xdata$Mes)))
for( i in levels(xdata$Mes)){
  #index variable
  j <-match(i, levels(xdata$Mes))
  nu <- nrow(xdata[xdata["Mes"] == i,])
  d <- rbind(xdata[xdata["Mes"] == i,], xdata[xdata["Mes"] == i,], xdata[xdata["Mes"] == i,])
  #creating one column that has the presence-absence of at for the gjam models and for the univariate
  #model the variable takes on another factor level (2)
  d$mode <- factor(c(rep(0, nu), rep(1, nu), rep(2, nu)), labels = c("GJAM with absent Anopheles", 
                                                                     "GJAM with present Anopheles",
                                                                     "Univariate Model"))
  #adding the predictions
  d$pred <- c(dr_gj_con$pred_gj_con[dr_gj_con$Mes == i], ggd$univariate[ggd$Mes == i])
  
  #Doing the ggplot
  plot <- ggplot(data = d, aes(x = NDVIBEF_2000, y = pred, shape = mode)) +
    geom_point() + ylab("Predicted Probability") + xlab("NDVIBEF_2000 in Standard Units") +
    ggtitle (paste0("Response Curve of NDVI Before for Culex perexiguus in GJAM \n Conditional on Presence of Anophles Atroparvus vs. Univariate Model for Month ", i)) +
    theme(legend.title = element_blank())
  plot
  plots_uv_con[[j]] <- plot
}
plots_uv_con[[6]]

#################################################################################################
##For Anopheles troparvus
###Response Curve für IA_500 :predict the probability of presence for the different ia-values,
#holding the other covariates at their mean (0) or at a constant factor level (mes)

#a sequence of IA_500-values from min to max with 50 steps
ia <- seq(min(df$IA_500), max(df$IA_500), length.out = 50)
#names of covariates
nc <- c("IA_500", "NDVIBEF_2000", "Mes")
#the mean of the covariates is 0, bc I normalized them
#creating  x-data dataframe with all the avareges (which are zero)
data <- as.data.frame(matrix(0, ncol = length(nc), nrow = length(ia)))
names(data) <- nc  
#replace IA_500 0s with the sequence from above
data$IA_500 <-ia
#Next, we want to create a dataframe with "data" for all the different months (bc we want to have
#a response curve for every month separately)
#initialize empty data frame that will be filled in the loop over the months
xdata <- as.data.frame(matrix(, ncol = length(nc), nrow = 0))
names(xdata) <- nc  
#run a loop over all the months and append month after month (with the same 
#covariates) to the dataframe xdata
for (i in levels(df$Mes)){
  d<- data
  d$Mes <- i
  xdata <- rbind(xdata, d)
}
#Converting Mes to a factor variable
xdata$Mes <- factor(xdata$Mes, levels =c("Abril", "Mayo", "Junio", "Julio", "Agosto",
                                         "Septiembre"))
#Run the posterior simulations with xdata, we want the predictions on the "probability scale"
#That's why we use posterior_linpred with transform = T
d_res <- posterior_linpred(fit_fin_at, transform = T, newdata = xdata, seed = 333)
#getting the predictions on the "probability scale, by taking the mean per 
#observation/column
univariate <- colMeans(d_res)
#adding these predictions to data frame
ggd <- cbind(xdata, univariate)
#maiking the plot
response_at_uv <- ggplot(data = ggd, aes(x = IA_500, y = univariate, color = Mes)) +
  geom_point() + ylab("Predicted Probability") + xlab("IA_500 in Standard Units") +
  ggtitle ("Response Curve of Inundation Area for Anopheles troparvus in Univariate Model") +
  guides(color=guide_legend(title="Month"))

##Doing the same thing for gjam
#for unconditional predictions
#Converting the Mes-Variable to a bunch of dummies as demanded by gjam
dr_gj <- gjamDummy(xdata$Mes, xdata)
#rename certain columns to coefficient names of joint_fin output
names(dr_gj)[1:2] <- c("IA500", "NDVIBEF2000")
#add an intercept to the data frame
dr_gj$intercept <- rep(1 , nrow(xdata)) 
#add already calculated quadratic terms
dr_gj$"I(IA500^2)" <- (dr_gj$IA500)^2
dr_gj$"I(NDVIBEF2000^2)" <- 0
#add a superfluous line of 1s (gjamPredict doesnt run, if there is no
#variation in the covariate) >> we will delete this row after making the 
#predictions
dr_gj <- rbind(dr_gj,1)
#delete Mes variable (also needed for gjamPredict to run)
dr_gj$Mes <- NULL
#define the modelling settings for posterior simulations with dr_gj
newdata <- list(xdata = dr_gj, nsim = 4000)
#calculating the in-sample predictions (simulations)
sim <- gjamPredict(output = joint_fin, newdata = newdata)

#take the first 300 rows as predictions (deleting the last row); Predictions are the means of the
#ychains
pred_gj_un <- sim$sdList$yMu[1:300,2] 
dr_gj <- dr_gj[-nrow(dr_gj),]
#readding the Mes column
dr_gj$Mes <- ggd$Mes
#making the ggplot dataframe
ggd_gj_un <- cbind(dr_gj, pred_gj_un)
response_at_un <- ggplot(data = ggd_gj_un, aes(x = IA500, y = pred_gj_un, color = Mes)) +
  geom_point() + ylab("Predicted Probability") + xlab("IA_500 in Standard Units") +
  ggtitle ("Response Curve of Inundation Area for Anopheles troparvus in GJAM") +
  guides(color=guide_legend(title="Month"))

##Plot gjam and rstanarm response curves in the same plot

#prepare the dataframe to feed into ggplot (just add predictions of gjam to dataframe of rstanarm)

ggd$multivariate <- ggd_gj_un$pred_gj_un

#Do the ggplot 

response_at_uv_un <- ggplot(data = ggd, aes(x = IA_500, color = Mes)) +
  geom_point(aes(y = univariate, shape = "univariate")) + 
  geom_point(aes(y = multivariate, shape = "multivariate")) + ylab("Predicted Probability") + 
  xlab("IA_500 in Standard Units") +
  ggtitle ("Univariate vs. Unconditional Multivariate Response Curves of Inundation Area for Anopheles troparvus") +
  guides(color=guide_legend(title="Month"), shape = guide_legend(title="Model"))

#Why are the response curves different? Does that contradict our hypothesis that unconditional
#predictions dont differ from univariate predictions? Why do they cross? (heterogeneous effect (
#of month depends on value of IA_500)) Why is for september multivariate probability more positive,
#but for all the other months the univariate model

#DOing the same thing with conditional predictions of gjam

#replicate dr_gj: We want to store both conditional prediction in one data frame (so we run gjamPredict
#only once and obtain both conditional predictions
dr_gj_con <- rbind(dr_gj, dr_gj)


#delete Mes variable (also needed for gjamPredict to run)
dr_gj_con$Mes <- NULL
#add a superfluous line of 1s (gjamPredict doesnt run, if there is no
#variation in the covariate) >> we will delete this row after making the 
#predictions
dr_gj_con <- rbind(dr_gj_con,1)
#create the PA-values for the Culex, the species we condition on (we create an addional one 
#so that it has the same length as dr_gj_con)
ydata <- tibble(c(rep(0, 300), rep(1, 301)))
names(ydata) <- "Cxperpre"


#define the modelling settings for posterior simulations with the conditional species as present
newdata <- list(xdata = dr_gj_con, ydataCond = ydata, nsim = 4000)

#calculating the in-sample predictions (simulations)
sim <- gjamPredict(output = joint_fin, newdata = newdata)

#take the first 600 rows as predictions (deleting the last row); Predictions are the means of the
#ychains
pred_gj_con <- sim$sdList$yMu[1:600,2] 

#adding the conditional presence of Anotrophes to our data storage data frame
dr_gj_con$cond <- factor(ydata$Cxperpre, labels = c("absent", "present"))
#deleting the last row of our data (remember: we added only bc this way the predict function runs)
dr_gj_con <- dr_gj_con[-601,]
#readding the Mes column
dr_gj_con$Mes <- ggd$Mes

dr_gj_con <- cbind(dr_gj_con, pred_gj_con)
response_at_gj_con <- ggplot(data = dr_gj_con, aes(x = IA500, y = pred_gj_con, color = Mes, shape = cond)) +
  geom_point() + ylab("Predicted Probability") + xlab("IA_500 in Standard Units") +
  ggtitle ("Response Curve of Inundation Area for Anopheles troparvus in GJAM \n Conditional on Presence of Anophles Atroparvus") +
  guides(color=guide_legend(title="Month"), shape = guide_legend(title = "Culex Perexigus"))

#looks pretty consistent (response is larger, if Culex is present)

###Doing a plot of gjam conditional predictions and univariate prediction (We do month by month)

#Preparing the input data
plots_uv_con <- vector(mode = "list", length = length(levels(xdata$Mes)))
for( i in levels(xdata$Mes)){
  #index variable
  j <-match(i, levels(xdata$Mes))
  nu <- nrow(xdata[xdata["Mes"] == i,])
  d <- rbind(xdata[xdata["Mes"] == i,], xdata[xdata["Mes"] == i,], xdata[xdata["Mes"] == i,])
  #creating one column that has the presence-absence of at for the gjam models and for the univariate
  #model the variable takes on another factor level (2)
  d$mode <- factor(c(rep(0, nu), rep(1, nu), rep(2, nu)), labels = c("GJAM with absent Culex", 
                                                                     "GJAM with present Culex",
                                                                     "Univariate Model"))
  #adding the predictions
  d$pred <- c(dr_gj_con$pred_gj_con[dr_gj_con$Mes == i], ggd$univariate[ggd$Mes == i])
  
  #Doing the ggplot
  plot <- ggplot(data = d, aes(x = IA_500, y = pred, shape = mode)) +
    geom_point() + ylab("Predicted Probability") + xlab("IA_500 in Standard Units") +
    ggtitle (paste0("Response Curve of Inundation Area for Anopheles troparvus in GJAM \n Conditional on Presence of Culex Perexigus vs. Univariate Model for Month ", i)) +
    theme(legend.title = element_blank())
  plots_uv_con[[j]] <- plot
}
plots_uv_con[[6]]

################################################
#for NDVIBEF#

#a sequence of NDVIBEF-values from min to max with 50 steps
nd <- seq(min(df$NDVIBEF_2000), max(df$NDVIBEF_2000), length.out = 50)
#names of covariates
nc <- c("IA_500", "NDVIBEF_2000", "Mes")
#the mean of the covariates is 0, bc I normalized them
#creating  x-data dataframe with all the avareges (which are zero)
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
for (i in levels(df$Mes)){
  d<- data
  d$Mes <- i
  xdata <- rbind(xdata, d)
}
#Converting Mes to a factor variable
xdata$Mes <- factor(xdata$Mes, levels =c("Abril", "Mayo", "Junio", "Julio", "Agosto",
                                         "Septiembre"))
#Run the posterior simulations with xdata, we want the predictions on the "probability scale"
#That's why we use posterior_linpred with transform = T
d_res <- posterior_linpred(fit_fin_at, transform = T, newdata = xdata, seed = 333)
#getting the predictions on the "probability scale, by taking the mean per 
#observation/column
univariate <- colMeans(d_res)
#adding these predictions to data frame
ggd <- cbind(xdata, univariate)
#maiking the plot
response_at_uv <- ggplot(data = ggd, aes(x = NDVIBEF_2000, y = univariate, color = Mes)) +
  geom_point() + ylab("Predicted Probability") + xlab("NDVIBEF_2000 in Standard Units") +
  ggtitle ("Response Curve of NDVI Before for Anopheles troparvus in Univariate Model") +
  guides(color=guide_legend(title="Month"))

#Interesting overall shape >> kinda contradicts the niche concept >> we have maxima at the edges
##Doing the same thing for gjam

#for unconditional predictions
#Converting the Mes-Variable to a bunch of dummies as demanded by gjam and storing whole data in
#dr_gj
dr_gj <- gjamDummy(xdata$Mes, xdata)
#rename certain columns to coefficient names of joint_fin output
names(dr_gj)[1:2] <- c("IA500", "NDVIBEF2000")
#add an intercept to the data frame
dr_gj$intercept <- rep(1 , nrow(xdata)) 
#add quadratic terms
dr_gj$"I(IA500^2)" <- 0
dr_gj$"I(NDVIBEF2000^2)" <- (dr_gj$NDVIBEF2000)^2
#add a superfluous line of 1s (gjamPredict doesnt run, if there is no
#variation in the covariate) >> we will delete this row after making the 
#predictions
dr_gj <- rbind(dr_gj,1)
#delete Mes variable (also needed for gjamPredict to run)
dr_gj$Mes <- NULL
#define the modelling settings for posterior simulations with dr_gj
newdata <- list(xdata = dr_gj, nsim = 4000)
#calculating the in-sample predictions (simulations)
sim <- gjamPredict(output = joint_fin, newdata = newdata)

#take the first 300 rows as predictions (deleting the last row); Predictions are the means of the
#ychains
pred_gj_un <- sim$sdList$yMu[1:300,2] 
dr_gj <- dr_gj[-nrow(dr_gj),]
#readding the Mes column
dr_gj$Mes <- ggd$Mes
#making the ggplot dataframe
ggd_gj_un <- cbind(dr_gj, pred_gj_un)
response_at_un <- ggplot(data = ggd_gj_un, aes(x = NDVIBEF2000, y = pred_gj_un, color = Mes)) +
  geom_point() + ylab("Predicted Probability") + xlab("NDVIBEF_2000 in Standard Units") +
  ggtitle ("Response Curve of NDVI Before for Anopheles troparvus in GJAM") +
  guides(color=guide_legend(title="Month"))

##Plot gjam and rstanarm response curves in the same plot

#prepare the dataframe to feed into ggplot (just add predictions of gjam to dataframe of rstanarm)

ggd$multivariate <- ggd_gj_un$pred_gj_un

#Do the ggplot 

response_at_uv_un <- ggplot(data = ggd, aes(x = NDVIBEF_2000, color = Mes)) +
  geom_point(aes(y = univariate, shape = "univariate")) + 
  geom_point(aes(y = multivariate, shape = "multivariate")) + ylab("Predicted Probability") + 
  xlab("NDVIBEF_2000 in Standard Units") +
  ggtitle ("Univariate vs. Multivariate Response Curves of NDVI Before for Anopheles troparvus") +
  guides(color=guide_legend(title="Month"), shape = guide_legend(title="Model"))

#They have a whole different shape >> I dont't think they should behave differently, but apparently
#they do o_o



#DOing the same thing with conditional predictions of gjam

#replicate dr_gj: We want to store both conditional prediction in one data frame (so we run gjamPredict
#only once and obtain both conditional predictions
dr_gj_con <- rbind(dr_gj, dr_gj)


#delete Mes variable (also needed for gjamPredict to run)
dr_gj_con$Mes <- NULL
#add a superfluous line of 1s (gjamPredict doesnt run, if there is no
#variation in the covariate) >> we will delete this row after making the 
#predictions
dr_gj_con <- rbind(dr_gj_con,1)
#create the PA-values for the Anopheles, the species we condition on (we create an addional one 
#so that it has the same length as dr_gj_con)
ydata <- tibble(c(rep(0, 300), rep(1, 301)))
names(ydata) <- "Cxperpre"


#define the modelling settings for posterior simulations with the conditional species as present
newdata <- list(xdata = dr_gj_con, ydataCond = ydata, nsim = 4000)

#calculating the in-sample predictions (simulations)
sim <- gjamPredict(output = joint_fin, newdata = newdata)

#take the first 600 rows as predictions (deleting the last row); Predictions are the means of the
#ychains
pred_gj_con <- sim$sdList$yMu[1:600,2] 

#adding the conditional presence of Culex to our data storage data frame
dr_gj_con$cond <- factor(ydata$Cxperpre, labels = c("absent", "present"))
#deleting the last row of our data (remember: we added only bc this way the predict function runs)
dr_gj_con <- dr_gj_con[-601,]
#readding the Mes column
dr_gj_con$Mes <- ggd$Mes

dr_gj_con <- cbind(dr_gj_con, pred_gj_con)
response_at_gj_con <- ggplot(data = dr_gj_con, aes(x = NDVIBEF2000, y = pred_gj_con, color = Mes, shape = cond)) +
  geom_point() + ylab("Predicted Probability") + xlab("NDVIBEF_2000 in Standard Units") +
  ggtitle ("Response Curve of NDVI Before for Anopheles troparvus in GJAM \n Conditional on Presence of Anophles Atroparvus") +
  guides(color=guide_legend(title="Month"), shape = guide_legend(title = "Anopheles atroparvus"))

#Interesting to note: Effect of conditional species is different depending on the month (But I 
#think that makes sense)

###Doing a plot of gjam conditional predictions and univariate prediction (We just take any month)


#Preparing the input data
plots_uv_con <- vector(mode = "list", length = length(levels(xdata$Mes)))
for( i in levels(xdata$Mes)){
  #index variable
  j <-match(i, levels(xdata$Mes))
  nu <- nrow(xdata[xdata["Mes"] == i,])
  d <- rbind(xdata[xdata["Mes"] == i,], xdata[xdata["Mes"] == i,], xdata[xdata["Mes"] == i,])
  #creating one column that has the presence-absence of at for the gjam models and for the univariate
  #model the variable takes on another factor level (2)
  d$mode <- factor(c(rep(0, nu), rep(1, nu), rep(2, nu)), labels = c("GJAM with absent Culex", 
                                                                     "GJAM with present Culex",
                                                                     "Univariate Model"))
  #adding the predictions
  d$pred <- c(dr_gj_con$pred_gj_con[dr_gj_con$Mes == i], ggd$univariate[ggd$Mes == i])
  
  #Doing the ggplot
  plot <- ggplot(data = d, aes(x = NDVIBEF_2000, y = pred, shape = mode)) +
    geom_point() + ylab("Predicted Probability") + xlab("NDVIBEF_2000 in Standard Units") +
    ggtitle (paste0("Response Curve of NDVI Before for Anopheles troparvus in GJAM \n Conditional on Presence of Culex perexiguus vs. Univariate Model for Month ", i)) +
    theme(legend.title = element_blank())
  plot
  plots_uv_con[[j]] <- plot
}
plots_uv_con[[1]]
#same pattern every month

#############################6.1d Variable Importance with DALEX################################
###################################### For Culex perexiguus
##################for univariate model

#Define our covariates as xdata
xdata <- train[,c("Mes", "IA_500", "NDVIBEF_2000")]
# create custom predict function for rstanarm "posterior_predict" function
#We take the mean of the simulations as the prediction  on the probability-scale. Reason: Easier
#to retrieve these values in gjam
pred_uv <- function(model, newdata)  {
  return(posterior_predict(model, newdata) %>% apply(2, mean))
}

#create the explain object (core object of DALEX) which contains the data and
#the predict function
dal_cp <- explain(fit_fin_cp, xdata, y = train$Cxperpre, predict_function = pred_uv,
                  type = "classification", label = "Univariate Probit")


# calculate the permutation-based variable importance measure (Here: the difference in
#(1-AUC) between original data and permuted data per covariate)
set.seed(1980)
vi_cp <- model_parts(dal_cp, type = "difference", B = 50)

#plot the results
plot(vi_cp) + 
  labs(title = "Variable Importance over 50 Permuations", subtitle = "created for the univariate probit model of Culex perexiguus") 
#So, Mes is the  most important variable, followed by IA_500 and at last 
#NDVIBEF_2000. When permuting the Mes variable entries and then predicting
#our response, the resulting AUC is roughly .11 worse than for the 
#predictions without the permutations. 


##################For gjam unconditional predictions

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
  #return the mean of the y-chains for every observation as the predictions
  return(pre$sdList$yMu[,1])
}

#make the explainer
dal_cpgj <- explain(joint_fin, xdata, y = train$Cxperpre,
                    predict_function = pred_gj, type = "classification",
                    label = "Unconditional Multivariate Probit")

#permutation-based variable importance measure
set.seed(1980)
vi_cpgj <- model_parts(dal_cpgj, type = "difference", B = 50) 

#plot the results
plot(vi_cpgj) + labs(title = "Variable Importance over 50 Permutations", subtitle = "created for the multivariate probit model of Culex perexiguus")
#Most important variable is IA_500, then Mes and then NDVIBEF_2000. This is slightly different than
#for the univariate case, but is the difference significant?

################## Conditional multivariate predictions

#plotting the two: 
#adding the other prediction types
d_gg <- tibble(var = c(vi_cp$variable, vi_cpgj$variable), loss = c(vi_cp$dropout_loss, vi_cpgj$dropout_loss), type = c(vi_cp$label, vi_cpgj$label))
#converting var and type to a factor
d_gg$var <- factor(d_gg$var, levels = c("NDVIBEF_2000", "IA_500", "Mes", "_baseline_", "_full_model_"))
d_gg$type <- factor(d_gg$type, levels = c("Univariate Probit", "Unconditional Multivariate Probit"))

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


#################################### Conditional multivariate predictions



#xdata with PA of Anopheles. This way, we treat the PA-data like a variable for which we can
#calculate VI
xdata_con <- as_tibble(cbind(xdata, y_train[,2]))

#make the predict function, including PA-anopheles as a covariate to
#permute/manipulate
pred_con <- function(model, newdata){
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
#Anopheles is the most important variable, followed by IA_500 and Mes. So, the environmental 
#covariates have the same order as in the unconditional case. But their magnitude is way smaller.
#The reason is that the PA-data takes away some of the predictive power of the environmental 
#covariates.


###########################################################################
######################################For Anopheles
################## For Univariate Probit Model
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
#our response, the resulting AUC is roughly .2 worse than for the 
#predictions without the permutations. 

##################For unconditional predictions in gjam

#We need to slightly change our predict function (just the last line) compared to the Culex case
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
#Most important variable is Mes, then NDVI_BEF and then IA, similar results as for the univariate
#model

#Plotting VI of both models in one plot

# making a dataframe with all the important data: 1. The variable which is permuted, 2. the drop 
#out loss for every run and 3. the according model type labels
d_ggat <- tibble(var = c(vi_at$variable, vi_atgj$variable),
                 loss = c(vi_at$dropout_loss, vi_atgj$dropout_loss),
                 type = c(vi_at$label, vi_atgj$label))
#converting var and type to a factor
d_ggat$var <- factor(d_ggat$var, levels = c("IA_500", "NDVIBEF_2000", "Mes", "_full_model_", 
                                            "_baseline_"))
d_ggat$type <- factor(d_ggat$type, levels = c("Univariate Probit", "Unconditional Multivariate Probit"))

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

##################Doing it for conditional predictions
#Treat the conditioning as a covariate (inlcude it as one variable in your
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
#Culex isnt that important o_o; is that in accordance with the AUCs? AUC is just a little over .01
#better for conditional model >> So, I guess so!
#Mes variable does not lose much of its importance measured in dif. The other two variables lose 
#considerably.

##############################################################################################
##################################################6.2. Out-of-Sample##########################
###############a. Conditional Predictions of gjam vs. Unconditional Predictions of############
#gjam vs. predictions of univariate model

#We evaluate predictive performance on our test set with the AUC
#####For Culex perexiguus
#We make the three predictions ((i) univariate, (ii) unconditional multivariate and (iii)
#conditional multivariate)

###(i) univariate predictions on the observation scale
pred_cp_uv <- posterior_predict(fit_fin_cp, newdata = test, seed = 333)
#take the average of these draws per observation as an estimation of the 
# "expected" predicted y-value/predicted probability of occurence
pred_cp_uv <- colMeans((pred_cp_uv))
#AUC
auc_cp_uv <- auc(response = test$Cxperpre, predictor = pred_cp_uv)
auc_cp_uv
#a AUC of .76 >> is ok for our purposes

###Unconditional multivariate
#predictions
#preparing the test data for gjam
# getting rid of all the unused variables
test_gj <- test[,c("Mes", "IA_500", "NDVIBEF_2000")] 
#month dummies
test_gj <- gjamDummy(test_gj$Mes, test_gj)
#adding intercept
test_gj$intercept <- rep(1 , nrow(test_gj)) 
#adding the squared terms
test_gj$"I(IA500^2)" <- (test_gj$IA_500)^2
test_gj$"I(NDVIBEF_2000^2)" <- (test_gj$NDVIBEF_2000)^2
#define the modelling settings for unconditional predictions with test_gj
newdata <- list(xdata = test_gj, nsim = 4000)
#calculating the in-sample predictions (simulations)
pred_cp_gj_un <- gjamPredict(output = joint_fin, newdata = newdata)
#retrieving the predictions (the mean of the simulations)
p_cp_gj_un <- pred_cp_gj_un$sdList$yMu[,1] 
#AUC
auc_cp_mvun <- auc(response = test$Cxperpre, predictor = p_cp_gj_un)
auc_cp_mvun
#AUC roughly .76 >> the same as the univariate model >> corroborates our hypothesis

###Plot the univariate conditions against the unconditional multivariate ones
d_gg_uc_cp <- data.frame(cbind(pred_cp_uv, p_cp_gj_un, y_test$Cxperpre))
names(d_gg_uc_cp) <- c("cp_uv", "cp_mv_un", "cp")
#make a factor out cp. This comes in handy for the plotting
d_gg_uc_cp$cp <- factor(d_gg_uc_cp$cp, levels = c(0, 1), labels = c("absent", "present"))

uvvsmvun_cp <- ggplot(d_gg_uc_cp, aes(x=cp_uv, y=cp_mv_un, shape=factor(cp))) + geom_point() +
  ggtitle("Univariate vs. Unconditional Multivariate Predictions for Culex perexiguus") +
  xlab("Predictions from Univariate Probit ") + 
  ylab("Unconditional Predictions from Multivariate Probit") +
  labs(shape = "Observed PA of Culex perexiguus") + ylim(0, 1)
uvvsmvun_cp
#They look as if they are centered around the identity line >> the predictions are more or less 
#the same! There is a slight tendency, though: For low predictions multivariate predictions are
#higher than univariate, but for high predictions univariate predictions are higher than 
#multivariate ones. This weakly contradicts hypothesis 3?

###Conditional predictions
#Culex perexiguus conditioned on Anopheles troparvus
#storing input data in newdata, inclduing the PA of Anopheles
set.seed(333)
newdata <- list(xdata = test_gj, ydataCond = y_test[,2], nsim = 4000) # conditionally predict out-of-sample
#Doing the actual prediction
pred_cp_mvco      <- gjamPredict(output = joint_fin, newdata = newdata)
#retrieving the predictions (the mean of the simulations)
p_cp_mvco <- pred_cp_mvco$sdList$yMu[,1] 
#AUC
auc_cp_mvco <- auc(response = test$Cxperpre, predictor = p_cp_mvco)
auc_cp_mvco
#AUC roughly .81 >> So, there is a considerable improvement of .05 in AUC >> corroborates
#our hypothesis that PA data enhances predictions :)

#Plotting Univariate Predictions vs. Conditional Multivariate Predictions
d_gg_cp <- data.frame(cbind(pred_cp_uv, p_cp_mvco, pred_cp_mvco$prPresent[,2], y_test$Cxperpre))
names(d_gg_cp) <- c("cp_uv", "cp_mv", "at", "cp")
#Make the PA-variables to factors
d_gg_cp$cp <- factor(d_gg_cp$cp, levels = c(0, 1), labels = c("absent", "present"))
d_gg_cp$at <- factor(d_gg_cp$at, levels = c(0, 1), labels = c("absent", "present"))

provsgj_cp <- ggplot(d_gg_cp, aes(x=cp_uv, y=cp_mv, color=at, shape = cp)) + geom_point() +
  ggtitle("Univariate vs. Conditional Multivariate Predictions for Culex perexiguus") +
  xlab("Predictions from Univariate Probit ") + ylab("Conditional Predictions from Multivariate Probit") +
  labs( color = "PA of Anopheles troparvus", shape = "PA of Culex perexiguus") + ylim(0,1)
provsgj_cp
#You can see that there is a positive linear relationship between the predictions
#of the two models grouped by the PA of Anopheles troparvus (the species we conditioned on).
#This indicates
#that both models roughly do the same/environmental signals are treated similarily.
#The slopes of the two lines are not equal to 1 (Would we expect this? I think, we do, if we 
#assume that the environmental coefficients are the same PUH, I still dont know the answer
#to this question). You can see that the conditioning on 
#Anapheles troparvus has a clear effect (The blue and red points form two distinct groups)
#on the multivariate predictions compared to the univariate predictions. The multivariate model
#predicts a roughly 12.5 % points
#higher probability of occurence for plots where Anopheles troparvus is present compared to
#plots where it's absent.



#####for Anopheles troparvus
#We evaluate predictive performance on our test set with the AUC
#We make the three predictions ((i) univariate, (ii) unconditional multivariate and (iii)
#conditional multivariate)

###(i) univariate predictions on the observation scale
pred_at_uv <- posterior_predict(fit_fin_at, newdata = test, seed = 333)
#take the mean of these draws per observation as an estimation of the 
# "expected" predicted y-value/predicted probability of occurence
pred_at_uv <- colMeans((pred_at_uv))
#AUC
auc_at_uv <- auc(response = test$Anatropre, predictor = pred_at_uv)
auc_at_uv
#a AUC of .85 >> is ok for our purposes

###Unconditional multivariate predictions

#define the modelling settings for unconditional predictions with test_gj
newdata <- list(xdata = test_gj, nsim = 4000)
#calculating the in-sample predictions (simulations)
pred_at_gj_un <- gjamPredict(output = joint_fin, newdata = newdata)
#retrieving the predictions (the mean of the simulations)
p_at_gj_un <- pred_at_gj_un$sdList$yMu[,2] 
#AUC
auc_at_mvun <- auc(response = test$Anatropre, predictor = p_at_gj_un)
auc_at_mvun
#AUC roughly .84 >> roughly the same as the univariate model >> corroborates our hypothesis

###Plot the univariate conditions against the unconditional multivariate ones
d_gg_uc_at <- data.frame(cbind(pred_at_uv, p_at_gj_un, y_test$Anatropre))
names(d_gg_uc_at) <- c("at_uv", "at_mv_un", "at")
d_gg_uc_at$at <- factor(d_gg_uc_at$at, levels = c(0, 1), labels = c("absent", "present"))

uvvsmvun_at <- ggplot(d_gg_uc_at, aes(x=at_uv, y=at_mv_un, shape=factor(at))) + geom_point() +
  ggtitle("Univariate vs. Unconditional Multivariate Predictions for Culex perexiguus") +
  xlab("Predictions from Univariate Probit ") + 
  ylab("Unconditional Predictions from Multivariate Probit") +
  labs(shape = "True PA of Culex perexiguus") + ylim(0, 1)
uvvsmvun_at
#Predictions are on one line, but not the identity line >> contradicts our hypothesis that 
#univariate probit predictions and unconditional multivariate probit predictions do not differ
#As for Culex, mv predictions are higher for low predictions and the other way around for high
#predictions

###Conditional predictions
#Anopheles troparvus conditioned on Culex perexiguus 
#storing input data in newdata, inclduing the PA of Culex
set.seed(333)
newdata <- list(xdata = test_gj, ydataCond = y_test[,1], nsim = 4000) # conditionally predict out-of-sample
#Doing the actual prediction
pred_at_mvco      <- gjamPredict(output = joint_fin, newdata = newdata)
#retrieving the predictions (the mean of the simulations)
p_at_mvco <- pred_at_mvco$sdList$yMu[,2] 
#AUC
auc_at_mvco <- auc(response = test$Anatropre, predictor = p_at_mvco)
auc_at_mvco
#AUC roughly .89 >> So, there is a considerable improvement of .04 in AUC >> corroborates
#our hypothesis that PA data enhances predictions :)

#Plotting Univariate Predictions vs. Conditional Multivariate Predictions
d_gg_at <- data.frame(cbind(pred_at_uv, p_at_mvco, y_test$Anatropre, y_test$Cxperpre))
names(d_gg_at) <- c("at_uv", "at_mv", "at", "cp")
#Make the PA-variables to factors
d_gg_at$cp <- factor(d_gg_at$cp, levels = c(0, 1), labels = c("absent", "present"))
d_gg_at$at <- factor(d_gg_at$at, levels = c(0, 1), labels = c("absent", "present"))

provsgj_at <- ggplot(d_gg_at, aes(x=at_uv, y=at_mv, color=cp, shape = at)) + geom_point() +
  ggtitle("Univariate vs. Conditional Multivariate Predictions for Anopheles troparvus") +
  xlab("Predictions from Univariate Probit ") + ylab("Conditional Predictions from Multivariate Probit") +
  labs( color = "PA of Culex perexiguus", shape = "PA of Anopheles troparvus") + ylim(0,1)
provsgj_at

#We can see that the conditioning on the PA of CUlex has a clear effect on predictions of 
#mv-probit. It's roughly .125 high, if Culex is present. Moreover, we detect the same pattern
#as in the unconditional vs. uv-plot: for low predictions mv has higher predictions and for 
#high predictions the other way around.

####################################6.2b Uncertainty in the predictions#######################
#We take a look at the dispersion of the predictions for the three prediction types.  

######For Culex perexiguus
###For the univariate model
pred_cpuv <- posterior_linpred(fit_fin_cp, transform = T, newdata = test, seed = 333)

#mean of the (standard deviations of the  predictions per observation)
sd_cpuv <- apply(pred_cpuv, 2, sd) %>% mean
#sd roughly .07
#for the unconditional predictions
sd_cp_mvun <- apply(pred_cp_gj_un$ychains[,1:140], 2, sd) %>% mean
#sd (.26) is a lot higher than for the univariate case
#for the conditional predictions
sd_cp_mvco <- apply(pred_cp_mvco$ychains[,1:140], 2, sd) %>% mean
#DAMN, the sd (.25) is way higher than in the univariate case >> would mean that the condtional 
#predictions are a lot more uncertain; slightly lower than the unconditional predictions.
#So, the uncertainty is probably induced by the different modelling approaches (multi vs. 
#univariate + difference in fitting between gjam and rstanarm), not the conditioning. I think,
#the reason is that gjam additonally (to the linear predictor with the environmental 
#covariates)draws the residual error from the multivariate normal distribution>> 
#More dispersion

#Plotting uncertainties of predictions with boxplots for single randomly drawn observations
i <- sample(seq_len(nrow(test)), size = 1)

#make the dataframe with all the predictions of observation i for all three prediction types
d_gg_bpso <- data.frame(c(pred_cpuv[,i], pred_cp_gj_un$ychains[,i], pred_cp_mvco$ychains[,i]))
names(d_gg_bpso) <- c("pred")
#adding a column with the prediction types
d_gg_bpso$type <- c(rep("Univariate", length(pred_cpuv[,i])),
                    rep("Unconditional Multivariate", length(pred_cp_gj_un$ychains[,i])),
                    rep("Conditional Multivariate", length(pred_cp_mvco$ychains[,i])))
d_gg_bpso$obs <- i
# Boxplot shows Boxes: 25% and 75 % Quantile; vertical line: median;
#lower whisker = smallest observation greater than or equal to lower hinge - 1.5 * IQR; 
#points: "remaining" outliers (die Erläuterungen komm in den Text unter der Abbildung)
ggplot(d_gg_bpso, aes(x = pred, y = type)) +
  geom_boxplot(aes(color = type)) + labs(title = paste0("Boxplot of 4000-simulated Predictions for Observation ", i),
                                         y = "Prediction Type",
                                         subtitle = "Culex perexiguus",
                                         x = "Prediction on the Probability Scale") + 
  theme(legend.position = "none") 
#>> We see that the univariate boxplot is much narrower

####95%-Confidence bands of predictions along predicted probabilities of all three predictions
#types

#calculate the .025 and .975 quantiles for each prediction type and each observation on the
#simulated predictions
quant_cpuv <- colQuantiles(pred_cpuv, probs=c(.025, .975))
quant_cpun <- colQuantiles(pred_cp_gj_un$ychains[,1:140], probs=c(.025, .975))
quant_cpco <- colQuantiles(pred_cp_mvco$ychains[,1:140], probs=c(.025, .975))

#make the data frame for the "cedribility intervall"confidence band" (correct word choice?)
#(ribbon argument of ggplot)
frame_cp <- data.frame(lower = c(quant_cpuv[,1], quant_cpun[,1], quant_cpco[,1]),
                       upper = c(quant_cpuv[,2], quant_cpun[,2], quant_cpco[,2]),
                       single = c(pred_cp_uv, p_cp_gj_un, p_cp_mvco),
                       type = factor(c(rep("Univariate Predictions", length(pred_cp_uv)),
                                       rep("Unconditional Multivariate Predictions", length(p_cp_gj_un)),
                                       rep("Conditional Multivariate Predictions", length(p_cp_mvco))),
                                     levels = c("Univariate Predictions", "Unconditional Multivariate Predictions",
                                                "Conditional Multivariate Predictions")))

#make the plot
ggplot(frame_cp, aes(single, single))+
  geom_line() +
  geom_ribbon(data=frame_cp,aes(ymin=lower,ymax=upper),alpha=0.3) +
  facet_wrap(~type, nrow =3)+ 
  labs(title = "95%-Confidence Bands of Predictions", subtitle = "Culex perexiguus",
       y = "Prediction", x = "Prediction")

#confidence bands for the the multivariate predictions are much, much broader (Hypothesis 2).
#Unconditional and condtional predictions seem to behave similarly

####################For Anopheles
###For the univariate model
pred_atuv <- posterior_linpred(fit_fin_at, transform = T, newdata = test, seed = 333)

#mean of the (standard deviations of the  predictions per observation)
sd_atuv <- apply(pred_atuv, 2, sd) %>% mean
#sd roughly .06
#for the unconditional predictions
sd_at_mvun <- apply(pred_at_gj_un$ychains[,141:280], 2, sd) %>% mean
#sd (.23) is a lot higher than for the univariate case
#for the conditional predictions
sd_at_mvco <- apply(pred_at_mvco$ychains[,141:280], 2, sd) %>% mean
#the sd (.22) is way higher than in the univariate case >> would mean that the condtional 
#predictions are a lot more uncertain; slightly lower than the unconditional predictions.
#So, the uncertainty is probably induced by the different modelling approaches (multi vs. 
#univariate + difference in fitting between gjam and rstanarm), not the conditioning. I think,
#the reason is that gjam additonally (to the linear predictor with the environmental 
#covariates)draws the residual error from the multivariate normal distribution>> 
#More dispersion

#Plotting uncertainties of predictions with boxplots for single randomly drawn observations

#make the dataframe with all the predictions of observation i for all three prediction types
d_gg_bpso <- data.frame(c(pred_atuv[,i], pred_at_gj_un$ychains[,140+i], pred_at_mvco$ychains[,140+i]))
names(d_gg_bpso) <- c("pred")
#adding a column with the prediction types
d_gg_bpso$type <- c(rep("Univariate", length(pred_atuv[,i])),
                    rep("Unconditional Multivariate", length(pred_at_gj_un$ychains[,i])),
                    rep("Conditional Multivariate", length(pred_at_mvco$ychains[,i])))
d_gg_bpso$obs <- i
# Boxplot shows Boxes: 25% and 75 % Quantile; vertical line: median;
#lower whisker = smallest observation greater than or equal to lower hinge - 1.5 * IQR; 
#points: "remaining" outliers (die Erläuterungen komm in den Text unter der Abbildung)
ggplot(d_gg_bpso, aes(x = pred, y = type)) +
  geom_boxplot(aes(color = type)) + labs(title = paste0("Boxplot of 4000-simulated Predictions for Observation ", i),
                                         subtitle = "Anopheles troparvus",
                                         y = "Prediction Type",
                                         x = "Prediction on the Probability Scale") + 
  theme(legend.position = "none") 
#>> We see that the univariate boxplot is much narrower

####95%-Confidence bands of predictions along predicted probabilities of all three predictions
#types

#calculate the .025 and .975 quantiles for each prediction type and each observation on the
#simulated predictions
quant_atuv <- colQuantiles(pred_atuv, probs=c(.025, .975))
quant_atun <- colQuantiles(pred_at_gj_un$ychains[,141:280], probs=c(.025, .975))
quant_atco <- colQuantiles(pred_at_mvco$ychains[,141:280], probs=c(.025, .975))

#make the data frame for the "cedribility intervall"confidence band" (correct word choice?)
#(ribbon argument of ggplot)
frame_at <- data.frame(lower = c(quant_atuv[,1], quant_atun[,1], quant_atco[,1]),
                       upper = c(quant_atuv[,2], quant_atun[,2], quant_atco[,2]),
                       single = c(pred_at_uv, p_at_gj_un, p_at_mvco),
                       type = factor(c(rep("Univariate Predictions", length(pred_at_uv)),
                                       rep("Unconditional Multivariate Predictions", length(p_at_gj_un)),
                                       rep("Conditional Multivariate Predictions", length(p_at_mvco))),
                                     levels = c("Univariate Predictions", "Unconditional Multivariate Predictions",
                                                "Conditional Multivariate Predictions")))

#make the plot
ggplot(frame_at, aes(single, single))+
  geom_line() +
  geom_ribbon(data=frame_at,aes(ymin=lower,ymax=upper),alpha=0.3) +
  facet_wrap(~type, nrow =3)+ 
  labs(title = "95%-Confidence Bands of Predictions", subtitle = "Anopheles troparvus",
       y = "Prediction", x = "Prediction")

#confidence bands for the the multivariate predictions are much, much broader (Hypothesis 2).
#Unconditional and condtional predictions seem to behave similarly,
