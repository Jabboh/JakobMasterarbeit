#I want to compare coefficients of the univariate and the multivariate models, viz. their size
#and their credibility intervals. I also want to create plots that visualize this.

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
library(bayesplot)
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

fit_cp <- stan_glm(Cxperpre ~ IA_500 + NDVI_500 + NDVIBEF_2000 + Mes, data = train, family = binomial(link = "probit"), init_r = 1.5, seed = 333)
summary(fit_cp)

# This is the plot provided by the rstanarm package >> looks good, but how do I put it together with
#the gjam plots?
plot(fit_cp) + 
  ggplot2::ggtitle("Posterior medians \n with 80% and 95% credible intervals")
plot(fit, pars = "size", regex_pars = "period", 
     ci_level = 0.95, outer_level = 1, show_density = TRUE)


###########GJAM Model
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
ml   <- list(ng = 4000, burnin = 100, typeNames = types)

#is 2000 for number of Gibbs steps ok?

#runnig GJAM
joint <- gjam(~ IA_500 + NDVI_500 + NDVIBEF_2000 + Mayo + Junio + Julio + Agosto + Septiembre, ydata = y_train, xdata = train_gj, modelList = ml)
summary(joint)
settings <- list(PLOTY = F, PLOTX = F, PREDICTX = F)
gjamPlot(joint, settings)
joint$parameters$betaMu

#This might do the trick: https://github.com/stan-dev/bayesplot/issues/232
#The posteriors of both models
posterior_1 <- as.matrix(fit_cp)
posterior_2 <- as.matrix(joint$chains$bgibbs[,1:9])
colnames(posterior_2) <- colnames(posterior_1)
combined <- rbind(mcmc_intervals_data(posterior_1), mcmc_intervals_data(posterior_2))
combined$model <- rep(c("Univariate Model", "Multivariate Model"), each = ncol(posterior_1))
theme_set(bayesplot::theme_default())
pos <- position_nudge(y = ifelse(combined$model == "Multivariate Model", 0, 0.1))
ggplot(combined, aes(x = m, y = parameter, color = model)) + 
  geom_point(position = pos) +
  geom_linerange(aes(xmin = ll, xmax = hh), position = pos)
