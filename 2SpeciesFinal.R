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

#For CP
waic_cp <- waic(fit_cp)

#Dropping IA_500
fit_cp1_ia <- stan_glm(Cxperpre ~ (NDVI_500 + IABEF_2000 + NDVIBEF_2000)^2 +
                     Mes + I(NDVI_500^2) + I(IABEF_2000^2) +
                     I(NDVIBEF_2000^2), data = train,
                   family = binomial(link = "probit"),init_r = .7, seed = 333)

waic_cp1_ia <- waic(fit_cp1_ia)

fit_cp1_ndvi <- stan_glm(Cxperpre ~ (IA_500 + IABEF_2000 + NDVIBEF_2000)^2 +
                     Mes + I(IA_500^2) + I(IABEF_2000^2) +
                     I(NDVIBEF_2000^2), data = train,
                   family = binomial(link = "probit"),init_r = .7, seed = 333)
waic_cp1_ndvi <- waic(fit_cp1_ndvi)

fit_cp1_iabef <- stan_glm(Cxperpre ~ (IA_500 + NDVI_500 + NDVIBEF_2000)^2 +
                     Mes + I(IA_500^2) + I(NDVI_500^2) +
                     I(NDVIBEF_2000^2), data = train,
                   family = binomial(link = "probit"),init_r = .7, seed = 333)
waic_cp1_iabef <- waic(fit_cp1_iabef)

fit_cp1_ndvibef <- stan_glm(Cxperpre ~ (IA_500 + NDVI_500 + IABEF_2000)^2 +
                     Mes + I(IA_500^2) + I(NDVI_500^2) + I(IABEF_2000^2),
                     data = train, family = binomial(link = "probit"),init_r = .7, 
                     seed = 333)
waic_cp1_ndvibef <- waic(fit_cp1_ndvibef)

fit_cp1_mes <- stan_glm(Cxperpre ~ (IA_500 + NDVI_500 + IABEF_2000 + NDVIBEF_2000)^2 +
                    I(IA_500^2) + I(NDVI_500^2) + I(IABEF_2000^2) +
                     I(NDVIBEF_2000^2), data = train,
                   family = binomial(link = "probit"),init_r = .7, seed = 333)
waic_cp1_mes <- waic(fit_cp1_mes)


loo_compare(waic_cp1_ia, waic_cp1_ndvi, waic_cp1_iabef, waic_cp1_ndvibef,
            waic_cp1_mes, waic_cp)
#The model with the highest elpd_diff is best model >> It's the model without the ndvi_500
#variable terms. Hence, we remove the ndvi. 

#Round 2
#Next, we work with the model without the month dummies and repeat the procedure
waic_cp1 <- waic(fit_cp1_ndvi)

#Dropping IA_500
fit_cp2_ia <- stan_glm(Cxperpre ~ (IABEF_2000 + NDVIBEF_2000)^2 +
                         Mes + I(IABEF_2000^2) +
                         I(NDVIBEF_2000^2), data = train,
                       family = binomial(link = "probit"),init_r = .7, seed = 333)

waic_cp2_ia <- waic(fit_cp2_ia)

fit_cp2_mes <- stan_glm(Cxperpre ~ (IA_500 + IABEF_2000 + NDVIBEF_2000)^2 +
                          I(IA_500^2) + I(IABEF_2000^2) +
                          I(NDVIBEF_2000^2), data = train,
                        family = binomial(link = "probit"),init_r = .7, seed = 333)
waic_cp2_mes <- waic(fit_cp2_mes)

fit_cp2_iabef <- stan_glm(Cxperpre ~ (IA_500 + NDVIBEF_2000)^2 +
                            Mes + I(IA_500^2) +
                            I(NDVIBEF_2000^2), data = train,
                          family = binomial(link = "probit"),init_r = .7, seed = 333)
waic_cp2_iabef <- waic(fit_cp2_iabef)

fit_cp2_ndvibef <- stan_glm(Cxperpre ~ (IA_500 + IABEF_2000)^2 +
                              Mes + I(IA_500^2) + I(IABEF_2000^2), data = train,
                            family = binomial(link = "probit"),init_r = .7, seed = 333)
waic_cp2_ndvibef <- waic(fit_cp2_ndvibef)

#comparing the WAICs
loo_compare(waic_cp2_ia, waic_cp2_mes, waic_cp2_iabef, waic_cp2_ndvibef,
            waic_cp1)
#The model with the highest elpd_diff is best model >> It's the model without the IABEF-terms
#terms. Hence, we remove the variable IABEF_2000 and all its associations. 

#Round 3
#Next, we work with the model without the ndvi_500 and the IABEF_2000 variable
#and repeat the procedure
waic_cp2 <- waic(fit_cp2_iabef)


fit_cp3_mes <- stan_glm(Cxperpre ~ (IA_500 + NDVIBEF_2000)^2 +
                           I(IA_500^2) +
                           I(NDVIBEF_2000^2), data = train,
                         family = binomial(link = "probit"),init_r = .7, seed = 333)
waic_cp3_mes <- waic(fit_cp3_mes)

fit_cp3_ia <- stan_glm(Cxperpre ~ (NDVIBEF_2000)^2 +
                            Mes +
                            I(NDVIBEF_2000^2), data = train,
                          family = binomial(link = "probit"),init_r = .7, seed = 333)
waic_cp3_ia <- waic(fit_cp3_ia)

fit_cp3_ndvibef <- stan_glm(Cxperpre ~ (IA_500)^2 +
                              Mes + I(IA_500^2), data = train,
                            family = binomial(link = "probit"),init_r = .7, seed = 333)
waic_cp3_ndvibef <- waic(fit_cp3_ndvibef)

#comparing the WAICs
loo_compare(waic_cp3_mes, waic_cp3_ia, waic_cp3_ndvibef,
            waic_cp2)
#The best model is the model with all the remaining covariates >> We do not remove 
# a covariate and all its associations


######this needs to be corrected!
#round 4 (1 interaction round)
#Next we check wether the interaction terms improve the WAIC. We always drop all the
#interaction terms associated with one variable

fit_cp4_ia <- stan_glm(Cxperpre ~ IA_500 + NDVIBEF_2000 +
                         Mes + I(IA_500^2) +
                         I(NDVIBEF_2000^2), data = train, refresh = 0,
                       family = binomial(link = "probit"),init_r = .7, seed = 333)
waic_cp4_ia <- waic(fit_cp4_ia)

fit_cp4_ndvibef <- stan_glm(Cxperpre ~ IA_500 + NDVIBEF_2000 +
                            Mes + I(IA_500^2) +
                            I(NDVIBEF_2000^2), data = train, refresh = 0,
                          family = binomial(link = "probit"),init_r = .7, seed = 333)

waic_cp4_ndvibef <- waic(fit_cp4_ndvibef)

#comparing the WAICs
loo_compare(waic_cp4_ia, waic_cp4_ndvibef,
            waic_cp2)
#The best model is the one with no  interactions, hence we stick to the
#model without the interactions.

#round 5 (first quadratic round)
#Next, we check whether the individual quadratic effects improve the WAIC. Hence, we 
#drop the quadratic terms one by one and decide as we did above

#our new baseline model



fit_cp5_ia <- stan_glm(Cxperpre ~ IA_500 + NDVIBEF_2000 +
                           Mes +
                           I(NDVIBEF_2000^2), data = train, refresh = 0,
                         family = binomial(link = "probit"),init_r = .7, seed = 333)
waic_cp5_ia <- waic(fit_cp5_ia)


fit_cp5_ndvibef <- stan_glm(Cxperpre ~ IA_500 + NDVIBEF_2000 +
                              Mes + I(IA_500^2), data = train, refresh = 0,
                            family = binomial(link = "probit"),init_r = .7, seed = 333)
waic_cp5_ndvibef <- waic(fit_cp5_ndvibef)


#comparing the WAICs
loo_compare(waic_cp5_ia, waic_cp5_ndvibef,waic_cp4_ia)
#The best model is the one with all the quadratic terms. Hence, we stick to the 
#model with all quadratic terms


#In conclusion, the best model (the most parsimonous with the best prediction
#performance for cp is the model with the linear specification of IA_500, NDVIBEF_2000,
# and mes as well as the quadratic terms of IA_500 and NDVIBEF_2000

fit_cp_b <- stan_glm(Cxperpre ~ IA_500 + NDVIBEF_2000 + Mes + I(IA_500^2) +
                           I(NDVIBEF_2000^2), data = train,
                         family = binomial(link = "probit"),init_r = .7, seed = 333)
