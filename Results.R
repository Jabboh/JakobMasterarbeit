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
#install.packages("devtools")

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
  ggtitle("Univariate vs. Multivariate Coefficients and their 95 % - Credibility Intervals \n for Anopheles Troparvus") + 
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
  ggtitle ("Univariate vs. Unconditional Multivariate Response Curves \n of Inundation Area for Culex perexiguus") +
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
  ggtitle ("Response Curve of Inundation Area for Anopheles Troparvus in Univariate Model") +
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
  ggtitle ("Response Curve of Inundation Area for Anopheles Troparvus in GJAM") +
  guides(color=guide_legend(title="Month"))

##Plot gjam and rstanarm response curves in the same plot

#prepare the dataframe to feed into ggplot (just add predictions of gjam to dataframe of rstanarm)

ggd$multivariate <- ggd_gj_un$pred_gj_un

#Do the ggplot 

response_at_uv_un <- ggplot(data = ggd, aes(x = IA_500, color = Mes)) +
  geom_point(aes(y = univariate, shape = "univariate")) + 
  geom_point(aes(y = multivariate, shape = "multivariate")) + ylab("Predicted Probability") + 
  xlab("IA_500 in Standard Units") +
  ggtitle ("Univariate vs. Multivariate Response Curves of Inundation Area for Anopheles Troparvus") +
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
  ggtitle ("Response Curve of Inundation Area for Anopheles Troparvus in GJAM \n Conditional on Presence of Anophles Atroparvus") +
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
    ggtitle (paste0("Response Curve of Inundation Area for Anopheles Troparvus in GJAM \n Conditional on Presence of Culex Perexigus vs. Univariate Model for Month ", i)) +
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
  ggtitle ("Response Curve of NDVI Before for Anopheles Troparvus in Univariate Model") +
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
  ggtitle ("Response Curve of NDVI Before for Anopheles Troparvus in GJAM") +
  guides(color=guide_legend(title="Month"))

##Plot gjam and rstanarm response curves in the same plot

#prepare the dataframe to feed into ggplot (just add predictions of gjam to dataframe of rstanarm)

ggd$multivariate <- ggd_gj_un$pred_gj_un

#Do the ggplot 

response_at_uv_un <- ggplot(data = ggd, aes(x = NDVIBEF_2000, color = Mes)) +
  geom_point(aes(y = univariate, shape = "univariate")) + 
  geom_point(aes(y = multivariate, shape = "multivariate")) + ylab("Predicted Probability") + 
  xlab("NDVIBEF_2000 in Standard Units") +
  ggtitle ("Univariate vs. Unconditional Multivariate Response Curves of NDVI Before for Anopheles Troparvus") +
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
  ggtitle ("Response Curve of NDVI Before for Anopheles Troparvus in GJAM \n Conditional on Presence of Anophles Atroparvus") +
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
    ggtitle (paste0("Response Curve of NDVI Before for Anopheles Troparvus in GJAM \n Conditional on Presence of Culex perexiguus vs. Univariate Model for Month ", i)) +
    theme(legend.title = element_blank())
  plot
  plots_uv_con[[j]] <- plot
}
plots_uv_con[[1]]
#same pattern every month