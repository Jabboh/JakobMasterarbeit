rm(list=ls())
setwd("C:\\Users\\jakob\\Documents\\JakobMasterarbeit\\Data")
#loading packages
library(readxl) #to read in the data
library(rstanarm) #Doing Bayesian probit GLMs for single species
library(pROC) #calculating the AUC
library(gjam) #Doing the joint estimation with the GJAM modelling aproach
library(ggplot2) #for plotting
library(DHARMa) # for checking in-sample validity of the models
library(dplyr)
####Data Preperation
#read in the data (Monthly species PA data for seven different mosquito species
#and according environmental covariates)
df <- read_excel("MonthlyData.xlsx")
str(df)
summary(df$Fecha) # Dates in 2006 dont make any sense. I assume that they put by accident 
#2006 instead of 2010. So I change these dates
df$Fecha <- as.POSIXct(sub("2006", "2010", df$Fecha))

#Reading in the spatial coordinates of the different trap locations
coords <- read_excel("Traps_coordenadas_geograficas.xls")

#Trap and Area mean the same thing, so we change the name in df from "area" to "trap"
names(df)[names(df)=="Area"] <- "trap"
#Canada is spelt differently, so I change the spelling in df to "Cañada"
df[,"trap"] <- lapply(df[,"trap"], gsub, pattern = "Ca?da", replacement = "Cañada", fixed = TRUE)

#adding lon-lat column to the data frame df
df <- merge(df, coords[, c("trap", "Norte", "Oeste")], by="trap", all.x= T, sort = F)

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

#Transform "Mes" (month when record was taken) into a factor variable
df$Mes <- factor(df$Mes, levels =c("Abril", "Mayo", "Junio", "Julio", "Agosto", "Septiembre"))
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

fit_cp <- stan_glm(Cxperpre ~ IA_500 + NDVI_500 + NDVIBEF_2000 + I(NDVI_500^2) + I(NDVIBEF_2000^2) + Mes, data = train, family = binomial(link = "probit"), init_r = 1.5, seed = 333)
summary(fit_cp)

#fitting the model for Anopheles atroparvus
fit_at <- stan_glm(Anatropre ~ IA_500 + NDVI_500 + NDVIBEF_2000 + I(NDVI_500^2) + I(NDVIBEF_2000^2) + Mes, data = train, family = binomial(link = "probit"), init_r = 1.5, seed = 333)
summary(fit_at)

###In-sample Validation


#Create a dharma object. For this, I specify the following arguments:
#1. 4000 simulations of fitted responses per observation
#(simulated response), 2. the observed responses (observedResponse),
# 3. the median of the 4000 simulated fitted responses which is the "expected" 
#value of the predicted y of an observation (fittedPredictedResponse),
#and other arguments concerning the specific handling of the "scaled residuals".
dharm_cp <- createDHARMa(simulatedResponse = t(posterior_predict(fit_cp)), observedResponse = fit_cp$y,
                         fittedPredictedResponse = posterior_predict(fit_cp) %>% apply(2, median), integerResponse = T, seed = 123,
                         method = "PIT")
#Why does posterior_predict product 0-1 values? Shouldn't it be the fitted y-values? 
plotQQunif(dharm_cp)
plot(dharm_cp)
#QQPLOT:
#Goal: Detect deviations from the expected distribution (here: uniform distribution)
#Here: looks like the observed distribution matches the expected distribution >> model
#correctly specified. Additionally, it gives us some p-values of statistical tests:
#1. KS test: test for correct distribution (Which distribution? of response?) >> here 
#not significant >> supports that I selected the correct distribution
#2. Dispersion test: no overdispersion
#3.Outliers test: no outliers (Note: values that are by default defined as values outside 
#the simulation envelope >> the more simulations, the fewer outliers)

#Boxplot:
#Scaled residuals vs. model predictions:
#We see tht residuals are for both posible predictions roughly the same

#Plot residuals against all covariates
plotResiduals(dharm_cp, train$Mes) #looks good >> no identication of quadratic effects
plotResiduals(dharm_cp, train$IA_500) #looks ok, no significant problems, but for .25 quantile
#no flat line >> could be indication of quadratic effect
plotResiduals(dharm_cp, train$NDVI_500)#looks periodic >> but also not 
#.75 quantile is still significant
plotResiduals(dharm_cp, train$NDVIBEF_2000)#looks ok
#looks good :), no longer quadratic
plotResiduals(dharm_cp, (train$NDVIBEF_2000)^2) # looks ok

#Ergo, quadratic terms might help???? Maybe compare it with the quadratic model?
#I mean the lines will never be super straight for all cases right?


hist(dharm_cp)
#looks pretty flat >> thumbs up! 

#different statistical tests
p <- testUniformity(dharm_cp)
p$p.value

# simulations, we have a lot >> no outliers
testDispersion(dharm_cp) # YOu should do this
#Calculating the residuals for both models
#guck dir alle Tests aus dem einen Teil in der Vignette nochmal an

##Test for temporal autocorrelation
dharm_cp_auto = recalculateResiduals(dharm_cp, group = train$Fecha)
plot(dharm_cp_auto)
testTemporalAutocorrelation(dharm_cp_auto, time =  unique(train$Fecha))
#does  seem to be a problem

#Spatial Autocorrelation

#I do not have the coordinates for trap "M29" --> one NA in the coordinates >> We need to remove that row
dharm_cp_spatial <- recalculateResiduals(dharm_cp, group = train$Norte)
testSpatialAutocorrelation(dharm_cp_spatial, 
                           x =  aggregate(train$Oeste, list(train$trap), mean)$x, 
                           y = aggregate(train$Norte, list(train$trap), mean)$x)
#No spatial autocorrelation >> Yeah!

###########DHARMA for anopheles

dharm_at <- createDHARMa(simulatedResponse = t(posterior_predict(fit_at)), observedResponse = fit_at$y,
                         fittedPredictedResponse = posterior_predict(fit_at) %>% apply(2, median), integerResponse = T, seed = 123,
                         method = "PIT")
#Why does posterior_predict product 0-1 values? Shouldn't it be the fitted y-values? 
plotQQunif(dharm_at)
plot(dharm_at)
#We got a problem! Und wie kann es sein, dass der dispersion test signifikant ist?
#QQPLOT:
#Goal: Detect deviations from the expected distribution (here: uniform distribution)
#Here: looks like the observed distribution matches the expected distribution >> model
#correctly specified. Additionally, it gives us some p-values of statistical tests:
#1. KS test: test for correct distribution (Which distribution? of response?) >> here 
#not significant >> supports that I selected the correct distribution
#2. Dispersion test: no overdispersion
#3.Outliers test: no outliers (Note: values that are by default defined as values outside 
#the simulation envelope >> the more simulations, the fewer outliers)

#Boxplot:
#Scaled residuals vs. model predictions:
#We see tht residuals are for both posible predictions roughly the same

#Plot residuals against all covariates
plotResiduals(dharm_at, train$Mes) 
#looks ok
plotResiduals(dharm_at, train$IA_500) #two quantiles are significant
plotResiduals(dharm_at, train$NDVI_500)
# ns, but looks periodic
plotResiduals(dharm_at, train$NDVIBEF_2000)
#looks good
plotResiduals(dharm_at, (train$NDVIBEF_2000)^2)
#looks ok
#Ergo, quadratic terms might help???? Maybe compare it with the quadratic model?
#I mean the lines will never be super straight for all cases right?
#or interactions? or random effects

hist(dharm_at)
#looks pretty flat >> thumbs up! 

#different statistical tests
p <- testUniformity(dharm_at)
p$p.value

testDispersion(dharm_at) # weird result. Did I do something wrong?

##Test for temporal autocorrelation
dharm_at_auto = recalculateResiduals(dharm_at, group = train$Fecha)
plot(dharm_at_auto)
testTemporalAutocorrelation(dharm_at_auto, time =  unique(train$Fecha))
#doesnt  seem to be a problem

#Spatial Autocorrelation

#I do not have the coordinates for trap "M29" --> one NA in the coordinates >> We need to remove that row
dharm_at_spatial <- recalculateResiduals(dharm_at, group = train$Norte)
testSpatialAutocorrelation(dharm_at_spatial, 
                           x =  aggregate(train$Oeste, list(train$trap), mean)$x, 
                           y = aggregate(train$Norte, list(train$trap), mean)$x)
#No spatial autocorrelation >> Yeah!

##########################Question
#Why does posterior_predict product 0-1 values? Shouldn't it be the fitted y-values?
#It makes sense to use the predicted values on the respnse scale, bc to calculate the scaled 
#residuals, we use the observed y-value whcih is 0 or 1. If, we would use the fitted y-values
#either all simulations would be greater than the observed y (for y = 0) or all
#smaller for y = 1. Hence, the scaled residuals would only take two values (0 or 1) >>
#it would no longer be a uniform distribution over the entire interval

#Why it works, even in the 0-1 case. Extreme case:Assume that the fuzzy term adds fuzziness of
#around .5 in one direction + Assume our model just predicts 1, then
#for all the 0s in the model the scale residual is between 0  and for all the 1s, 
#the scaled residual is between 0 and 1. Thus, the distribution is no longer uniform
#but skewed towards 0 ... Is this correct?
#library(matlab)
#x <- zeros(325, 4000)
#x[1,] <- 1

#dharm_test <- createDHARMa(simulatedResponse = x, observedResponse = fit_cp$y,
#fittedPredictedResponse = posterior_predict(fit_cp) %>% apply(2, median), integerResponse = T, seed = 123,
#method = "PIT")
#For some fucked up reason, it does not work
#Wie wäre es mit Schlucken? Ich glaube es macht auch Sinn, sagt meine Intuition ;)
######################Fitting a joint model of both species with GJAM

#Define the model settings
types <- c("PA", "PA") #variable types of the responses
s <- length(types) # number of species
#y-data to train the model
y_train <- y[train_id, ]

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
#define model/algorithm parameters
ml   <- list(ng = 2000, burnin = 100, typeNames = types)

#is 2000 for number of Gibbs steps ok?

#runnig GJAM
joint <- gjam(~ IA_500 + NDVI_500 + NDVIBEF_2000 + Mayo + Junio + Julio + Agosto + Septiembre, ydata = y_train, xdata = train_gj, modelList = ml)
summary(joint)

############################DHARMa-Analysis
#Create a dharma object. For this, I specify the following arguments:
#1. 4000 simulations of fitted responses per observation
#(simulated response), 2. the observed responses (observedResponse),
# 3. the median of the 4000 simulated fitted responses which is the "expected" 
#value of the predicted y of an observation (fittedPredictedResponse),
#and other arguments concerning the specific handling of the "scaled residuals".

#Prepare a newdata object for the gjam-posterior predictions
#add intercept to train data (for some reason we need to this, otherwise gjamPredict
#won't accept it as xdata)
train_gj$intercept <- rep(1 , nrow(train)) #adding the intercept
# getting rid of all the unused variables
train_gj <- train_gj[,c("intercept", "IA_500", "NDVI_500", "NDVIBEF_2000", "Mayo", "Junio", "Julio",
                        "Agosto", "Septiembre", "Fecha", "trap", "Norte", "Oeste")] 
# getting rid of all the unused variables

newdata <- list(xdata = train_gj, nsim = 4000)
#calculating the in-sample predictions (simulations)
sim <- gjamPredict(output = joint, newdata = newdata)
#I am pretty sure that the single draws from the predictive posterior are stored in ychains
#I am, however, a little confused how gjamPredict calculates the one-value predictions
#(prPresent). It's pretty close to the median
#For the above reasons I have to convert the fitted y-values to the observation scale.
#I think, it should just be >.5 -> 1, 0 otherwise
sims <- ifelse(sim$ychains>.5, 1, 0)
dharm_gj_un <- createDHARMa(simulatedResponse = t(sims),
                            observedResponse = append(fit_cp$y, fit_at$y, after = length(fit_cp$y)),
                            fittedPredictedResponse = append(sim$prPresent[,1], sim$prPresent[,2],
                                                             after = length(sim$prPresent[,1])),
                            integerResponse = T, seed = 123,method = "PIT")

plot(dharm_gj_un)
#looks pretty good
#different statistical tests
p <- testUniformity(dharm_gj_un)
p$p.value


###Covariates
plotResiduals(dharm_gj_un, rep(train$Mes,2)) 
#looks pretty good
plotResiduals(dharm_gj_un, rep(train$IA_500,2)) #looks goood
plotResiduals(dharm_gj_un, rep(train$NDVI_500, 2))
#.5 and .75 quantile is significant
plotResiduals(dharm_gj_un, rep(train$NDVIBEF_2000, 2))

##Test for temporal autocorrelation, we need to append FEcha with Fecha, bc we have
#S*n observations
dharm_gj_un_auto <- recalculateResiduals(dharm_gj_un, group = rep(train_gj$Fecha, 2))
plot(dharm_gj_un_auto)
testTemporalAutocorrelation(dharm_gj_un_auto, time =  unique(train$Fecha))
#doesn't  seem to be a problem


#####Spatial Autocorrelation

#adding lon-lat column to the training data

dharm_gj_spatial <- recalculateResiduals(dharm_gj_un, group = rep(train$trap, 2))
testSpatialAutocorrelation(dharm_gj_spatial, 
                           x =  aggregate(rep(train_gj$Oeste, 2), list(rep(train_gj$trap, 2)), mean)$x, 
                           y = aggregate(rep(train_gj$Norte, 2),
                                         list(rep(train_gj$trap, 2)), mean)$x)


#############
