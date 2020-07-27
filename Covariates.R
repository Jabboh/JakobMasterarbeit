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
######Plotting Covariates against the response
#Since my response is categorical, scatterplots will not help us to eye-ball functional relationships
#betweeen our predictors and the response. I implement a moving average of the response with respect
#to the predictors (Presences/(Presences/Absences) with respect to covariate)

#The span-argument in loess defines the size of the neighborhood for which the average is calculated
#>> this can potentially have a large impact on the curve (shape of the curve is used to eye-ball
#functional relationship)
####CP

#IA_500
cp_ia <- loess(Cxperpre ~ IA_500,data=df)
plot(Cxperpre ~ IA_500,data=df ,pch=19,cex=0.1)
j <- order(df$IA_500)
lines(df$IA_500[j],cp_ia$fitted[j],col="red",lwd=3)
#looks like a negative relationship, looks quadratic, but one could argue that it's locally linear

#NDVI_500
cp_nd <- loess(Cxperpre ~ NDVI_500,data=df, span = .75)
plot(Cxperpre ~ NDVI_500,data=df,pch=19,cex=0.1)
j <- order(df$NDVI_500)
lines(df$NDVI_500[j],cp_nd$fitted[j],col="red",lwd=3)
#looks fucking crazy, higher order polynomial... BUT for larger spans (e.g. 1), relationship looks
#quadratic

#NDVI Before_2000
cp_nb <- loess(Cxperpre ~ NDVIBEF_2000,data=df, span = 1)
plot(Cxperpre ~ NDVIBEF_2000,data=df,pch=19,cex=0.1)
j <- order(df$NDVIBEF_2000)
lines(df$NDVIBEF_2000[j],cp_nb$fitted[j],col="red",lwd=3)
#looks quadratic


####Interaction Effects
int <- df$NDVIBEF_2000 * df$NDVI_500
cp_ndnb <- loess(df$Cxperpre ~ int)
plot(df$Cxperpre ~ int,pch=19,cex=0.1)
j <- order(int)
lines(int[j],cp_ndnb$fitted[j],col="red",lwd=3)
summary(cp_nb)

####Eyeballing: Interaction Terms
#https://stackoverflow.com/questions/41988812/plot-smoothed-average-of-third-variable-by-x-and-y
train %>% group_by(IA_500 = cut(IA_500, 7),
                   NDVI_500 = cut(NDVI_500, 7)) %>%
  summarise(y=mean(Cxperpre)) %>%
  ggplot(aes(IA_500, NDVI_500, fill=y)) +
  geom_tile() + 
  theme_classic()
