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