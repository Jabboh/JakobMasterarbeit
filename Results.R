#############################Results of my analyis#############################################
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
#install.packages("matrixStats")

#loading packages
library(readxl) #to read in the data
library(rstanarm) #Doing Bayesian probit GLMs for single species
library(pROC) #calculating the AUC
library(gjam) #Doing the joint estimation with the GJAM modelling aproach
library(ggplot2) #for plotting
library(DHARMa) # for checking in-sample validity of the models
library(dplyr) # for simplified syntax and neater code
library(loo) #to calculate WAIC
library(bayesplot) #Some handy features for plotting in the Bayesian realm
library(matrixStats)#for quantile calculations
library(ggplot2)#for plotting
library(DALEX)#For Variable Importance

##########################################1.Data Preperation#######################################
#read in the data (Monthly species PA data for seven different mosquito species
#and according environmental covariates)
df <- read_excel("MonthlyData.xlsx")

#Checking the rough structure of the data set
str(df)
summary(df$Fecha) # Dates in 2006 dont make any sense. I assume that they put by
#accident 2006 instead of 2010. So I change these dates.
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

#Selecting our two species (perexiguus & Anopheles troparvus) for the analysis.
#extract the PA data for all species
spec <- df[,7:14]
#deleting the An_atroparvus column, because this is not PA-data
spec[,"An_atroparvus"] <- NULL
#we select Cxperpre and Anatropre for our first analysis:
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

#Fitting final univariate models

#for Culex perexiguus

fit_fin_cp <- stan_glm(Cxperpre ~ Mes + IA_500 + NDVIBEF_2000 + I(IA_500^2) + I(NDVIBEF_2000^2),
                       data = train, family = binomial(link = "probit"), init_r = 1.4,  seed = 333)

#for anopheles atropovarus
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
##########################################6. Results#############################################
##########################################6.1. In-sample:########################################
#####################################6.1.a. Comparing the coefficients: Size, Significance and credibility intervals###########
#Plotting all the coefficients in one plot per species

####For Culex perexiguus
#https://github.com/stan-dev/bayesplot/issues/232
#The posteriors of both models
#rstanarm takes the median as the parameter estimate
#We also take the median in the following graphs
posterior_1 <- as.matrix(fit_fin_cp) 
posterior_2 <- as.matrix(joint_fin$chains$bgibbsUn[,1:10])
colnames(posterior_2) <- colnames(posterior_1)

#create an object that stores the median and the two boundaries for our credibility interval
combined <- rbind(mcmc_intervals_data(posterior_1, prob_outer = .95), mcmc_intervals_data(posterior_2, prob_outer = .95))
combined$model <- rep(c("Univariate Model", "Multivariate Model"), each = ncol(posterior_1))
#prob_outer defines the credibility interval in our plot (here .95)

#plot the parameters of both models and the respective credibility intervals
theme_set(bayesplot::theme_default())
pos <- position_nudge(y = ifelse(combined$model == "Multivariate Model", 0, 0.1))
ggplot(combined, aes(x = m, y = parameter, color = model)) + 
  geom_point(position = pos) +
  geom_vline(xintercept = 0, linetype="dotted", color = "black", size=.5) +
  geom_errorbar(aes(xmin = ll, xmax = hh), position = pos, width = .1) +
  ggtitle("Univariate vs. Multivariate Coefficients and their 95 % - Credibility Intervals \n for Culex perexiguus") + 
  xlab("Value") + ylab("Coefficient") + labs(color="Model") 
#They look pretty much the same 

####for Anopheles atropovarus
#posteriors of both models
posterior_1 <- as.matrix(fit_fin_at)
posterior_2 <- as.matrix(joint_fin$chains$bgibbsUn[,-(1:10)])
colnames(posterior_2) <- colnames(posterior_1)
#create an object with median and 95% credibility intervals of the posterior
combined <- rbind(mcmc_intervals_data(posterior_1, prob_outer = .95), mcmc_intervals_data(posterior_2, prob_outer = .95))
combined$model <- rep(c("Univariate Model", "Multivariate Model"), each = ncol(posterior_1))

#plotting parameters with credibility intervals
pos <- position_nudge(y = ifelse(combined$model == "Multivariate Model", 0, 0.1))
coef_at <- ggplot(combined, aes(x = m, y = parameter, color = model)) + 
  geom_point(position = pos) +
  geom_errorbar(aes(xmin = ll, xmax = hh), position = pos, width = .1) +
  geom_vline(xintercept = 0, linetype="dotted", color = "black", size=.5) +
  ggtitle("Univariate vs. Multivariate Coefficients and their 95 % - Credibility Intervals \n for Anopheles Troparvus") + 
  xlab("Value") + ylab("Coefficient") + labs(color="Model") 
#They look pretty much the same 

#set plotting them back to default
theme_set(theme_grey())

#### the hard numbers: Tables with the coefficients and SEs

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

#####################################6.1.b Correlation between the responses (raw vs. residual correlation)##########
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
#dummies are different)! But there is still a lot of unexplained Correlation.


#####################################6.1.c. Response Curves#########################################
####For Culex perexiguus

###Response Curve für IA_500:predict the probability of presence for different ia-values 
#(over the range of ia-values in the dataset), holding the other covariates at their mean (0)
#or at a constant factor level (mes)

#a sequence of IA_500-values from min to max with 50 steps
ia <- seq(min(df$IA_500), max(df$IA_500), length.out = 50)
#names of covariates
nc <- c("IA_500", "NDVIBEF_2000", "Mes")
#the mean of the covariates is 0, bc we normalized them
#creating  x-data dataframe with all the averages (which are zero)
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

#Run the posterior simulations of the univariate model with xdata, we want the predictions on
#the "probability scale". That's why we use posterior_linpred with transform = T
d_res <- posterior_predict(fit_fin_cp, newdata = xdata, seed = 333)
#taking the mean of the posterior simulations as our prediction
univariate <- colMeans(d_res)
#adding these predictions to data frame
ggd <- cbind(xdata, univariate)
#maiking the plot "Response Curve of Inundation Area for Culex perexiguus in Univariate Model"
ggplot(data = ggd, aes(x = IA_500, y = univariate, color = Mes)) +
  geom_point() + ylab("Predicted Probability") + xlab("IA_500 in Standard Units") +
  ggtitle ("Response Curve of Inundation Area for Culex perexiguus in Univariate Model") +
  guides(color=guide_legend(title="Month"))

####Doing the same thing for the multivariate model

####for unconditional predictions
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
pred_gj_un <- sim$sdList$yMu[1:300,1] 
dr_gj <- dr_gj[-nrow(dr_gj),]
#readding the Mes column
dr_gj$Mes <- ggd$Mes
#making the ggplot dataframe
ggd_gj_un <- cbind(dr_gj, pred_gj_un)
#The plot "Response Curve of Inundation Area for Culex perexiguus with Unconditional Predictions \n in Multivariate Model"
ggplot(data = ggd_gj_un, aes(x = IA500, y = pred_gj_un, color = Mes)) +
  geom_point() + ylab("Predicted Probability") + xlab("IA_500 in Standard Units") +
  ggtitle ("Response Curve of Inundation Area for Culex perexiguus with Unconditional Predictions \n in Multivariate Model") +
  guides(color=guide_legend(title="Month"))

#Plot uni- and multivariate response curves in the same plot

#prepare the dataframe to feed into ggplot (just add predictions of multivariate to dataframe of univariate)

ggd$multivariate <- ggd_gj_un$pred_gj_un

#Do the ggplot 

ggplot(data = ggd, aes(x = IA_500, color = Mes)) +
  geom_point(aes(y = univariate, shape = "univariate")) + 
  geom_point(aes(y = multivariate, shape = "multivariate")) + ylab("Predicted Probability") + 
  xlab("IA_500 in Standard Units") +
  ggtitle ("Univariate vs. Unconditional Multivariate Response Curves \n of Inundation Area for Culex perexiguus") +
  guides(color=guide_legend(title="Month"), shape = guide_legend(title="Model"))

#Why are the response curves different? Does that contradict our hypothesis that unconditional
#predictions dont differ from univariate predictions? Why do they cross? Increasing IA has
#a stronger negative effect on predicted probabilities in univariate models.  Why is for
#september the predicted prob of multivariate model more positive, but for all the other months
#the predicted probs of univariate models?

#####Doing the same thing with conditional predictions of gjam

#replicate dr_gj: We want to store both conditional prediction (the one where Anopheles is 
#absent and the one where it is present) in one data frame (so that we run gjamPredict
#only once and obtain both conditional predictions)
dr_gj_con <- rbind(dr_gj, dr_gj)


#delete Mes variable (also needed for gjamPredict to run)
dr_gj_con$Mes <- NULL
#add a superfluous line of 1s (gjamPredict doesnt run, if there is no
#variation in the covariate) >> we will delete this row after making the 
#predictions
dr_gj_con <- rbind(dr_gj_con,1)
#create the PA-values for the Anopheles, the species we condition on (we create an additional
#one so that it has the same length as dr_gj_con)
ydata <- tibble(c(rep(0, 300), rep(1, 301)))
names(ydata) <- "Anatropre"


#define the modelling settings for posterior simulations with the conditional species as ydataCond
newdata <- list(xdata = dr_gj_con, ydataCond = ydata, nsim = 4000)

#calculating the in-sample predictions (simulations)
sim <- gjamPredict(output = joint_fin, newdata = newdata)

#take the first 600 rows as predictions (deleting the last row); Predictions are the means of the
#ychains
pred_gj_con <- sim$sdList$yMu[1:600,1] 

#adding the conditional presence of Anopheles to our data storage data frame
dr_gj_con$cond <- factor(ydata$Anatropre, labels = c("absent", "present"))
#deleting the last row of our data (remember: we added only bc this way the predict function runs)
dr_gj_con <- dr_gj_con[-601,]
#readding the Mes column
dr_gj_con$Mes <- ggd$Mes

#combining covariates and conditional predictions
dr_gj_con <- cbind(dr_gj_con, pred_gj_con)
#The plot:"Response Curve of Inundation Area for Culex perexiguus in GJAM Conditional on Presence of Anophles troparvus"
ggplot(data = dr_gj_con, aes(x = IA500, y = pred_gj_con, color = Mes, shape = cond)) +
  geom_point() + ylab("Predicted Probability") + xlab("IA_500 in Standard Units") +
  ggtitle ("Response Curve of Inundation Area for Culex perexiguus in GJAM \n Conditional on Presence of Anophles troparvus") +
  guides(color=guide_legend(title="Month"), shape = guide_legend(title = "Anopheles troparvus"))

#strange results, ordering of the effect of month on pred, changes depending on whether AT is present
#or not. Does that make sense? Is that even possible? >> I actually do not think so (But only september
#is different!)

####Doing a plot of gjam conditional predictions and univariate prediction per month

#make a list of plots with one plot for every month. Each plot contains response three response
#curves: (i) univariate model, (ii) conditional predictions with absent Anopheles and (iii)
#conditional predictions with present Anopheles
#Preparing the input data
plots_uv_con <- vector(mode = "list", length = length(levels(xdata$Mes)))
for( i in levels(xdata$Mes)){
  #index variable
  j <-match(i, levels(xdata$Mes))
  #Number of observations
  nu <- nrow(xdata[xdata["Mes"] == i,])
  #three-times the xdata in one data frame, so that every prediction type is represented
  d <- rbind(xdata[xdata["Mes"] == i,], xdata[xdata["Mes"] == i,], xdata[xdata["Mes"] == i,])
  #creating one column that has the presence-absence of Anopheles for the gjam models and for the univariate
  #model the variable takes on another factor level (2)
  d$mode <- factor(c(rep(0, nu), rep(1, nu), rep(2, nu)), labels = c("GJAM with absent Anopheles", 
                                                                     "GJAM with present Anopheles",
                                                                     "Univariate Model"))
  #adding the predictions
  d$pred <- c(dr_gj_con$pred_gj_con[dr_gj_con$Mes == i], ggd$univariate[ggd$Mes == i])
  
  #Doing the ggplot
  plot <- ggplot(data = d, aes(x = IA_500, y = pred, shape = mode)) +
    geom_point() + ylab("Predicted Probability") + xlab("IA_500 in Standard Units") +
    ggtitle (paste0("Response Curve of Inundation Area for Culex perexiguus in GJAM \n Conditional on Presence of Anophles troparvus vs. Univariate Model for Month ", i)) +
    theme(legend.title = element_blank())
  plots_uv_con[[j]] <- plot
}
#Plot for Junio
plots_uv_con[[3]]


########for NDVIBEF

####for univariate predictions
#a sequence of NDVIBEF-values from min to max with 50 steps
nd <- seq(min(df$NDVIBEF_2000), max(df$NDVIBEF_2000), length.out = 50)

#the mean of the covariates is 0, bc I normalized them
#creating  x-data dataframe with all the averages (which are zero)
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
#Run the posterior simulations with xdata
d_res <- posterior_predict(fit_fin_cp, newdata = xdata, seed = 333)
#getting the predictions on the "probability scale, by taking the mean per 
#observation/column
univariate <- colMeans(d_res)
#adding these predictions to the data frame
ggd <- cbind(xdata, univariate)
#maiking the plot: "Response Curve of NDVI Before for Culex perexiguus in Univariate Model"
ggplot(data = ggd, aes(x = NDVIBEF_2000, y = univariate, color = Mes)) +
  geom_point() + ylab("Predicted Probability") + xlab("NDVIBEF_2000 in Standard Units") +
  ggtitle ("Response Curve of NDVI Before for Culex perexiguus in Univariate Model") +
  guides(color=guide_legend(title="Month"))

#Interesting overall shape >> kinda contradicts the niche concept >> we have maxima at the edges

####Doing the same thing for gjam

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

#Plot: "Response Curve of NDVI Before for Culex perexiguus in GJAM"
ggplot(data = ggd_gj_un, aes(x = NDVIBEF2000, y = pred_gj_un, color = Mes)) +
  geom_point() + ylab("Predicted Probability") + xlab("NDVIBEF_2000 in Standard Units") +
  ggtitle ("Response Curve of NDVI Before for Culex perexiguus in GJAM") +
  guides(color=guide_legend(title="Month"))

##Plot gjam and rstanarm response curves in the same plot

#prepare the dataframe to feed into ggplot (just add multivariate predictions to dataframe of 
#univariate predictions)

ggd$multivariate <- ggd_gj_un$pred_gj_un

#Plot: "Univariate vs. Multivariate Response Curves of NDVI Before for Culex perexiguus"

ggplot(data = ggd, aes(x = NDVIBEF_2000, color = Mes)) +
  geom_point(aes(y = univariate, shape = "univariate")) + 
  geom_point(aes(y = multivariate, shape = "multivariate")) + ylab("Predicted Probability") + 
  xlab("NDVIBEF_2000 in Standard Units") +
  ggtitle ("Univariate vs. Multivariate Response Curves of NDVI Before for Culex perexiguus") +
  guides(color=guide_legend(title="Month"), shape = guide_legend(title="Model"))

#They have a  different shape >> I dont't think they should behave differently, but apparently
#they do.



####Doing the same thing with conditional predictions of gjam

#replicate dr_gj: We want to store both conditional prediction (the ones with present Anophles
#and the one with absent Anopheles) in one data frame (so we run gjamPredict
#only once and obtain both conditional predictions)
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


#define the modelling settings for posterior simulations
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

#adding predictions to dataframe dr_gj_con (our xdata)
dr_gj_con <- cbind(dr_gj_con, pred_gj_con)
ggplot(data = dr_gj_con, aes(x = NDVIBEF2000, y = pred_gj_con, color = Mes, shape = cond)) +
  geom_point() + ylab("Predicted Probability") + xlab("NDVIBEF_2000 in Standard Units") +
  ggtitle ("Response Curve of NDVI Before for Culex perexiguus in GJAM \n Conditional on Presence of Anophles troparvus") +
  guides(color=guide_legend(title="Month"), shape = guide_legend(title = "Anopheles troparvus"))

#strange results, ordering of the effect of month on pred, changes depending on whether
#Anopheles is present or not. Does that make sense? Is that even possible? 
#>> I actually do not think so (BUt only septemberis different!)

####Doing a plot of gjam conditional predictions and univariate prediction per month

#make a list of plots with one plot for every month. Each plot contains response three response
#curves: (i) univariate model, (ii) conditional predictions with absent Anopheles and (iii)
#conditional predictions with present Anopheles

#Preparing the input data
plots_uv_con <- vector(mode = "list", length = length(levels(xdata$Mes)))
for( i in levels(xdata$Mes)){
  #index variable
  j <-match(i, levels(xdata$Mes))
  #number of observation per month
  nu <- nrow(xdata[xdata["Mes"] == i,])
  #stacking covariates 3-times (for every prediction type)
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
    ggtitle (paste0("Response Curve of NDVI Before for Culex perexiguus in GJAM \n Conditional on Presence of Anophles troparvus vs. Univariate Model for Month ", i)) +
    theme(legend.title = element_blank())
  plot
  plots_uv_con[[j]] <- plot
}

#Example: Septiembre
plots_uv_con[[6]]


####For Anopheles ####
###Response Curve für IA_500 :predict the probability of presence for the different ia-values,
#holding the other covariates at their mean (0) or at a constant factor level (mes)

#a sequence of IA_500-values from min to max with 50 steps
ia <- seq(min(df$IA_500), max(df$IA_500), length.out = 50)

#the mean of the covariates is 0, bc I normalized them
#creating  x-data dataframe with all the avareges (which are zero)
data <- as.data.frame(matrix(0, ncol = length(nc), nrow = length(ia)))
names(data) <- nc  
#replace IA_500 0s with the sequence from above
data$IA_500 <-ia
#Next, we want to create a dataframe with "data" for all the different months (because we want
#to have a response curve for every month separately)
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
#Run the posterior simulations with xdata
d_res <- posterior_predict(fit_fin_at, newdata = xdata, seed = 333)
#getting the predictions on the "probability scale, by taking the mean per 
#observation/column
univariate <- colMeans(d_res)
#adding these predictions to data frame
ggd <- cbind(xdata, univariate)
#maiking the plot "Response Curve of Inundation Area for Anopheles Troparvus in Univariate Model"
ggplot(data = ggd, aes(x = IA_500, y = univariate, color = Mes)) +
  geom_point() + ylab("Predicted Probability") + xlab("IA_500 in Standard Units") +
  ggtitle ("Response Curve of Inundation Area for Anopheles Troparvus in Univariate Model") +
  guides(color=guide_legend(title="Month"))

####Doing the same thing for multivariate model
####for unconditional predictions
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
#variation in the covariate) >> we will delete this row after making the predictions
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
#plot "Response Curve of Inundation Area for Anopheles Troparvus in GJAM"
ggplot(data = ggd_gj_un, aes(x = IA500, y = pred_gj_un, color = Mes)) +
  geom_point() + ylab("Predicted Probability") + xlab("IA_500 in Standard Units") +
  ggtitle ("Response Curve of Inundation Area for Anopheles Troparvus in GJAM") +
  guides(color=guide_legend(title="Month"))

####Plot multi- and univariate response curves in the same plot

#prepare the dataframe to feed into ggplot (just add multivariate predictions to dataframe of 
#univariate predictions)
ggd$multivariate <- ggd_gj_un$pred_gj_un

#Plot "Univariate vs. Multivariate Response Curves of Inundation Area for Anopheles Troparvus"
ggplot(data = ggd, aes(x = IA_500, color = Mes)) +
  geom_point(aes(y = univariate, shape = "univariate")) + 
  geom_point(aes(y = multivariate, shape = "multivariate")) + ylab("Predicted Probability") + 
  xlab("IA_500 in Standard Units") +
  ggtitle ("Univariate vs. Multivariate Response Curves of Inundation Area for Anopheles Troparvus") +
  guides(color=guide_legend(title="Month"), shape = guide_legend(title="Model"))

#Why are the response curves different? Does that contradict our hypothesis that unconditional
#predictions dont differ from univariate predictions? Why do they cross? 

####Doing the same thing with conditional predictions of multivariate model

#replicate dr_gj: We want to store both conditional predictions (the one with present Culex
#and the one with absent Culex) in one data frame
dr_gj_con <- rbind(dr_gj, dr_gj)


#delete Mes variable (also needed for gjamPredict to run)
dr_gj_con$Mes <- NULL
#add a superfluous line of 1s (gjamPredict doesnt run, if there is no
#variation in the covariate) >> we will delete this row after making the predictions
dr_gj_con <- rbind(dr_gj_con,1)
#create the PA-values for  Culex, the species we condition on (we create an addional one 
#so that it has the same length as dr_gj_con)
ydata <- tibble(c(rep(0, 300), rep(1, 301)))
names(ydata) <- "Cxperpre"


#define the modelling settings for posterior simulations 
newdata <- list(xdata = dr_gj_con, ydataCond = ydata, nsim = 4000)

#calculating the in-sample predictions (simulations)
sim <- gjamPredict(output = joint_fin, newdata = newdata)

#take the first 600 rows as predictions (deleting the last row); Predictions are the means of the
#ychains
pred_gj_con <- sim$sdList$yMu[1:600,2] 

#adding the PA of Culex to our data storage data frame
dr_gj_con$cond <- factor(ydata$Cxperpre, labels = c("absent", "present"))
#deleting the last row of our data (remember: we added only bc this way the predict function 
#runs)
dr_gj_con <- dr_gj_con[-601,]
#readding the Mes column
dr_gj_con$Mes <- ggd$Mes

#making data frame for ggplot
dr_gj_con <- cbind(dr_gj_con, pred_gj_con)
#Plotting "Response Curve of Inundation Area for Anopheles Troparvus in GJAM Conditional on Presence of Anophles troparvus"
ggplot(data = dr_gj_con, aes(x = IA500, y = pred_gj_con, color = Mes, shape = cond)) +
  geom_point() + ylab("Predicted Probability") + xlab("IA_500 in Standard Units") +
  ggtitle ("Response Curve of Inundation Area for Anopheles Troparvus in GJAM \n Conditional on Presence of Anophles troparvus") +
  guides(color=guide_legend(title="Month"), shape = guide_legend(title = "Culex perexiguus"))

#looks pretty consistent (response is larger, if Culex is present)

####Doing a plot of gjam conditional predictions and univariate prediction (We do month by month)

#make a list of plots with one plot for every month. Each plot contains response three response
#curves: (i) univariate model, (ii) conditional predictions with absent Culex and (iii)
#conditional predictions with present Culex
#Preparing the input data
plots_uv_con <- vector(mode = "list", length = length(levels(xdata$Mes)))
for( i in levels(xdata$Mes)){
  #index variable
  j <-match(i, levels(xdata$Mes))
  #number of observations per month
  nu <- nrow(xdata[xdata["Mes"] == i,])
  #appending the covariates
  d <- rbind(xdata[xdata["Mes"] == i,], xdata[xdata["Mes"] == i,], xdata[xdata["Mes"] == i,])
  #creating one column that has the PA of Culex for the multivariate models and for the 
  #univariate model the variable takes on another factor level (2)
  d$mode <- factor(c(rep(0, nu), rep(1, nu), rep(2, nu)), labels = c("GJAM with absent Culex", 
                                                                     "GJAM with present Culex",
                                                                     "Univariate Model"))
  #adding the predictions
  d$pred <- c(dr_gj_con$pred_gj_con[dr_gj_con$Mes == i], ggd$univariate[ggd$Mes == i])
  
  #Doing the ggplot
  plot <- ggplot(data = d, aes(x = IA_500, y = pred, shape = mode)) +
    geom_point() + ylab("Predicted Probability") + xlab("IA_500 in Standard Units") +
    ggtitle (paste0("Response Curve of Inundation Area for Anopheles Troparvus in GJAM \n Conditional on Presence of Culex perexiguus vs. Univariate Model for Month ", i)) +
    theme(legend.title = element_blank())
  plots_uv_con[[j]] <- plot
}

#Example: July
plots_uv_con[[4]]
#Why do the curves have different forms?
####for NDVIBEF

###For the univariate model
#a sequence of NDVIBEF-values from min to max with 50 steps
nd <- seq(min(df$NDVIBEF_2000), max(df$NDVIBEF_2000), length.out = 50)
#names of covariates
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
#Run the posterior simulations with xdata
d_res <- posterior_predict(fit_fin_at, newdata = xdata, seed = 333)
#getting the predictions on the "probability scale, by taking the mean per 
#observation/column
univariate <- colMeans(d_res)
#adding these predictions to data frame
ggd <- cbind(xdata, univariate)
#making the plot "Response Curve of NDVI Before for Anopheles Troparvus in Univariate Model"
ggplot(data = ggd, aes(x = NDVIBEF_2000, y = univariate, color = Mes)) +
  geom_point() + ylab("Predicted Probability") + xlab("NDVIBEF_2000 in Standard Units") +
  ggtitle ("Response Curve of NDVI Before for Anopheles Troparvus in Univariate Model") +
  guides(color=guide_legend(title="Month"))

#Interesting overall shape >> kinda contradicts the niche concept >> we have maxima at the edges

####Doing the same thing for multivariate model
####for unconditional predictions

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

#take the first 300 rows as predictions (deleting the last row) of second row (corresponding to
#Anopheles); Predictions are the means of the ychains
pred_gj_un <- sim$sdList$yMu[1:300,2] 
dr_gj <- dr_gj[-nrow(dr_gj),]
#readding the Mes column
dr_gj$Mes <- ggd$Mes
#making the ggplot dataframe
ggd_gj_un <- cbind(dr_gj, pred_gj_un)
#Plot "Response Curve of NDVI Before for Anopheles Troparvus in GJAM"
ggplot(data = ggd_gj_un, aes(x = NDVIBEF2000, y = pred_gj_un, color = Mes)) +
  geom_point() + ylab("Predicted Probability") + xlab("NDVIBEF_2000 in Standard Units") +
  ggtitle ("Response Curve of NDVI Before for Anopheles Troparvus in GJAM") +
  guides(color=guide_legend(title="Month"))

####Plot multi- and univariate response curves in the same plot

#prepare the dataframe to feed into ggplot (just add predictions of multivariate model
#to dataframe of univariate predictions)
ggd$multivariate <- ggd_gj_un$pred_gj_un

#Plot "Univariate vs. Unconditional Multivariate Response Curves of NDVI Before for Anopheles Troparvus" 
ggplot(data = ggd, aes(x = NDVIBEF_2000, color = Mes)) +
  geom_point(aes(y = univariate, shape = "univariate")) + 
  geom_point(aes(y = multivariate, shape = "multivariate")) + ylab("Predicted Probability") + 
  xlab("NDVIBEF_2000 in Standard Units") +
  ggtitle ("Univariate vs. Unconditional Multivariate Response Curves of NDVI Before for Anopheles Troparvus") +
  guides(color=guide_legend(title="Month"), shape = guide_legend(title="Model"))

#They have a whole different shape >> I dont't think they should behave differently, but apparently
#they do.

#Doing the same thing with conditional predictions of gjam

#replicate dr_gj: We want to store both conditional predictions (the ones with present Culex and 
#the one with absent Culex) in one data frame 
dr_gj_con <- rbind(dr_gj, dr_gj)
#delete Mes variable (also needed for gjamPredict to run)
dr_gj_con$Mes <- NULL
#add a superfluous line of 1s (gjamPredict doesnt run, if there is no
#variation in the covariate) >> we will delete this row after making the predictions
dr_gj_con <- rbind(dr_gj_con,1)
#create the PA-values for Culex, the species we condition on (we create an addional row 
#so that it has the same length as dr_gj_con)
ydata <- tibble(c(rep(0, 300), rep(1, 301)))
names(ydata) <- "Cxperpre"

#define the modelling settings for posterior simulations with conditional species
newdata <- list(xdata = dr_gj_con, ydataCond = ydata, nsim = 4000)

#calculating the in-sample predictions (simulations)
sim <- gjamPredict(output = joint_fin, newdata = newdata)

#take the first 600 rows as predictions (deleting the last row); Predictions are the means of 
#the ychains
pred_gj_con <- sim$sdList$yMu[1:600,2] 

#adding the conditional presence of Culex to our data storage data frame
dr_gj_con$cond <- factor(ydata$Cxperpre, labels = c("absent", "present"))
#deleting the last row of our data (remember: we added only bc this way the predict function runs)
dr_gj_con <- dr_gj_con[-601,]
#readding the Mes column
dr_gj_con$Mes <- ggd$Mes

#creating ggplot dataframe
dr_gj_con <- cbind(dr_gj_con, pred_gj_con)

#Plot "Response Curve of NDVI Before for Anopheles Troparvus in GJAM Conditional on Presence of Anophles troparvus"
ggplot(data = dr_gj_con, aes(x = NDVIBEF2000, y = pred_gj_con, color = Mes, shape = cond)) +
  geom_point() + ylab("Predicted Probability") + xlab("NDVIBEF_2000 in Standard Units") +
  ggtitle ("Response Curve of NDVI Before for Anopheles Troparvus in GJAM \n Conditional on Presence of Anophles troparvus") +
  guides(color=guide_legend(title="Month"), shape = guide_legend(title = "Anopheles troparvus"))

#makes sense to me 

####Doing a plot of multivariate conditional predictions and univariate predictions (We do month by month)

#make a list of plots with one plot for every month. Each plot contains response three response
#curves: (i) univariate model, (ii) conditional predictions with absent Culex and (iii)
#conditional predictions with present Culex

#Preparing the input data
plots_uv_con <- vector(mode = "list", length = length(levels(xdata$Mes)))
for( i in levels(xdata$Mes)){
  #index variable
  j <-match(i, levels(xdata$Mes))
  #number of observations per month
  nu <- nrow(xdata[xdata["Mes"] == i,])
  d <- rbind(xdata[xdata["Mes"] == i,], xdata[xdata["Mes"] == i,], xdata[xdata["Mes"] == i,])
  #creating one column that has the presence-absence of Culex for the gjam models and for the univariate
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

#Example: April
plots_uv_con[[1]]
#####################################6.1d Variable Importance with DALEX################################
####For Culex perexiguus
####for univariate model

#Define our covariates as xdata
xdata <- train[,c("Mes", "IA_500", "NDVIBEF_2000")]
# create custom predict function for rstanarm "posterior_predict" function
#We take the mean of the simulations as the prediction  on the probability-scale. Reason: Easier
#to retrieve these values in gjam
pred_uv <- function(model, newdata)  {
  return(posterior_predict(model, newdata) %>% apply(2, mean))
}

#create the explain object (core object of DALEX) which contains the model, the data and
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

####For multivariate unconditional predictions

# create custom predict function, this is a little trickier, because we need
#to feed it xdata with the three columns (Mes, IA_500 and NDVIBEF) (bc we 
#want to get the variable importance of these three variables), but we cannot use gjamPredict
#with only these three xdata columns >> We need to manipulate the data first, before applying
#gjamPredict. We do this in our custumory gjam-predict function:
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


####plotting the two: 
#making a tibble with (i) a variable column indicating which variable is permuted, (ii) a loss column giving the 
#(1-Auc)-losses and (iii) a type column indicating the model type
d_gg <- tibble(var = c(vi_cp$variable, vi_cpgj$variable), loss = c(vi_cp$dropout_loss, vi_cpgj$dropout_loss), type = c(vi_cp$label, vi_cpgj$label))
#converting var and type to a factor
d_gg$var <- factor(d_gg$var, levels = c("NDVIBEF_2000", "IA_500", "Mes", "_baseline_", "_full_model_"))
d_gg$type <- factor(d_gg$type, levels = c("Univariate Probit", "Unconditional Multivariate Probit"))

#deleting baseline and full model (We do not need these loss values)
d_gg <-d_gg[!(d_gg$var == "_full_model_" | d_gg$var == "_baseline_"),]

#Plot "Variable Importance Boxplots"
ggplot(d_gg, aes(x = loss, y = var)) +
  geom_boxplot(aes(color = type)) + # Boxplot shows Boxes: 25% and 75 % Quantile; vertical line: median; lower whisker = smallest observation greater than or equal to lower hinge - 1.5 * IQR
  facet_wrap(~type, nrow =3) +
  labs(title = "Variable Importance Boxplots",
       subtitle = "Created for multi- and univariate models of Culex perexiguus with 50 permutation runs", y = "Variable",
       x = "(1 - AUC)-Loss after Permutations") +
  # Suppress the legend since color isn't actually providing any information
  theme(legend.position = "none") 

#### Conditional multivariate predictions

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
  #make a newdata list for gjam with ydataCond for the condtitioning species
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
#Doing it for 50 permutations
vi_cp_con <- model_parts(dal_cp_con, type = "difference", B = 5)

#plotting the results
plot(vi_cp_con) + labs(title = "Variable Importance", subtitle = "created for the multivariate probit model of Culex perexiguus conditional on Anopheles troparvus")

#Anopheles is the most important variable, followed by IA_500 and Mes. So, the environmental 
#covariates have the same order as in the unconditional case. But their magnitude is way smaller.
#The reason is that the PA-data takes away some of the predictive power of the environmental 
#covariates.

####For Anopheles####
###For Univariate Probit Model
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

####For unconditional predictions in gjam

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
  #return the mean of the y-chains for every observation, this time of at (column 2)
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

###Plotting VI of both models in one plot

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

#Plot the comparison
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
  #make a newdata list for gjam with ydataCond
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
#Doing it for 50 permutations
vi_at_con <- model_parts(dal_at_con, type = "difference", B = 50)

#plotting the results
plot(vi_at_con) + labs(title = "Variable Importance", subtitle = "created for the multivariate probit model of Anopheles troparvus conditional on Culex perexiguus")
#Culex isnt that important; is that in accordance with the AUCs? AUC is just a little over .01
#better for conditional model >> So, I guess so!
#Mes variable does not lose much of its importance measured in dif. The other two variables lose 
#considerably.

##########################################6.2. Out-of-Sample##########################
#####################################6.2.a. Conditional Predictions of gjam vs. Unconditional Predictions of gjam vs. predictions of univariate model####

#We evaluate predictive performance on our test set with the AUC.

#####For Culex perexiguus
#We make the three predictions ((i) univariate, (ii) unconditional multivariate and (iii)
#conditional multivariate)

###(i) univariate predictions on the observation scale
pred_cp_uv <- posterior_predict(fit_fin_cp, newdata = test, seed = 333)
#take the average of these draws per observation as an estimation of the 
# "expected" predicted y-value/predicted probability of occurence (prediction scale)
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
#make a factor out of cp. This comes in handy for the plotting
d_gg_uc_cp$cp <- factor(d_gg_uc_cp$cp, levels = c(0, 1), labels = c("absent", "present"))

#Plot "Univariate vs. Unconditional Multivariate Predictions for Culex perexiguus"
ggplot(d_gg_uc_cp, aes(x=cp_uv, y=cp_mv_un, shape=factor(cp))) + geom_point() +
  ggtitle("Univariate vs. Unconditional Multivariate Predictions for Culex perexiguus") +
  xlab("Predictions from Univariate Probit ") + 
  ylab("Unconditional Predictions from Multivariate Probit") +
  labs(shape = "Observed PA of Culex perexiguus") + ylim(0, 1)

#They look as if they are centered around the identity line >> the predictions are more or less 
#the same! There is a slight tendency, though: For low predictions multivariate predictions are
#higher than univariate, but for high predictions univariate predictions are higher than 
#multivariate ones. This weakly contradicts hypothesis 3?

###Conditional predictions
#Culex perexiguus conditioned on Anopheles troparvus
#storing input data in newdata, inclduing the PA of Anopheles
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

###Plotting Univariate Predictions vs. Conditional Multivariate Predictions
d_gg_cp <- data.frame(cbind(pred_cp_uv, p_cp_mvco, pred_cp_mvco$prPresent[,2], y_test$Cxperpre))
names(d_gg_cp) <- c("cp_uv", "cp_mv", "at", "cp")
#Make the PA-variables to factors
d_gg_cp$cp <- factor(d_gg_cp$cp, levels = c(0, 1), labels = c("absent", "present"))
d_gg_cp$at <- factor(d_gg_cp$at, levels = c(0, 1), labels = c("absent", "present"))

ggplot(d_gg_cp, aes(x=cp_uv, y=cp_mv, color=at, shape = cp)) + geom_point() +
  ggtitle("Univariate vs. Conditional Multivariate Predictions for Culex perexiguus") +
  xlab("Predictions from Univariate Probit ") + ylab("Conditional Predictions from Multivariate Probit") +
  labs( color = "PA of Anopheles troparvus", shape = "PA of Culex perexiguus") + ylim(0,1)

#You can see that there is a positive linear relationship between the predictions
#of the two models grouped by the PA of Anopheles troparvus (the species we conditioned on).
#This indicates that both models roughly do the same/environmental signals are treated 
#similarily. The slopes of the two lines are not equal to 1 (Would we expect this? I think,
#we do, if we assume that the environmental coefficients are the same PUH, I still dont know
#the answerto this question). You can see that the conditioning on Anopheles troparvus has a 
#clear effect (The blue and red points form two distinct groups) on the multivariate 
#predictions compared to the univariate predictions. The multivariate model predicts a roughly
#12.5 % points higher probability of occurence for plots where Anopheles troparvus is present 
#compared to plots where it's absent.

####for Anopheles ####
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
#Preparing the data
d_gg_uc_at <- data.frame(cbind(pred_at_uv, p_at_gj_un, y_test$Anatropre))
names(d_gg_uc_at) <- c("at_uv", "at_mv_un", "at")
d_gg_uc_at$at <- factor(d_gg_uc_at$at, levels = c(0, 1), labels = c("absent", "present"))
#Plotting "Univariate vs. Unconditional Multivariate Predictions for Anopheles troparvus"
ggplot(d_gg_uc_at, aes(x=at_uv, y=at_mv_un, shape=factor(at))) + geom_point() +
  ggtitle("Univariate vs. Unconditional Multivariate Predictions for Anopheles troparvus") +
  xlab("Predictions from Univariate Probit ") + 
  ylab("Unconditional Predictions from Multivariate Probit") +
  labs(shape = "True PA of Culex perexiguus") + ylim(0, 1)

#Predictions are on one line, but not the identity line >> contradicts our hypothesis that 
#univariate probit predictions and unconditional multivariate probit predictions do not differ
#Like for Culex, multivariate predictions are higher for low predictions and the other way
#around for high predictions.

###Conditional predictions
#Anopheles troparvus conditioned on Culex perexiguus 
#storing input data in newdata, inclduing the PA of Culex
newdata <- list(xdata = test_gj, ydataCond = y_test[,1], nsim = 4000) # conditionally predict out-of-sample
#Doing the actual prediction
pred_at_mvco      <- gjamPredict(output = joint_fin, newdata = newdata)
#retrieving the predictions (the mean of the simulations)
p_at_mvco <- pred_at_mvco$sdList$yMu[,2] 
#AUC
auc_at_mvco <- auc(response = test$Anatropre, predictor = p_at_mvco)
auc_at_mvco
#AUC roughly .89 >> So, there is a considerable improvement of .04 in AUC >> corroborates
#our hypothesis that PA-data of other species enhances predictions.

#Plotting Univariate Predictions vs. Conditional Multivariate Predictions
d_gg_at <- data.frame(cbind(pred_at_uv, p_at_mvco, y_test$Anatropre, y_test$Cxperpre))
names(d_gg_at) <- c("at_uv", "at_mv", "at", "cp")
#Make the PA-variables to factors
d_gg_at$cp <- factor(d_gg_at$cp, levels = c(0, 1), labels = c("absent", "present"))
d_gg_at$at <- factor(d_gg_at$at, levels = c(0, 1), labels = c("absent", "present"))

ggplot(d_gg_at, aes(x=at_uv, y=at_mv, color=cp, shape = at)) + geom_point() +
  ggtitle("Univariate vs. Conditional Multivariate Predictions for Anopheles troparvus") +
  xlab("Predictions from Univariate Probit ") + ylab("Conditional Predictions from Multivariate Probit") +
  labs( color = "PA of Culex perexiguus", shape = "PA of Anopheles troparvus") + ylim(0,1)


#We can see that the conditioning on the PA of CUlex has a clear effect on predictions of 
#mv-probit. It's roughly .125 high, if Culex is present. Moreover, we detect the same pattern
#as in the unconditional vs. uv-plot: for low predictions mv has higher predictions and for 
#high predictions the other way around.

#####################################6.2.b Uncertainty in the predictions#######################
#We take a look at the dispersion of the predictions for the three prediction types.  

####For Culex perexiguus
###For the univariate model; We need to take the predictions on the probability scale
#We obtain them, getting the posterior predictions of the linear predictor (posterior_linpred)
#and transforming these value with the link function (transform = T)
pred_cpuv <- posterior_linpred(fit_fin_cp, transform = T, newdata = test, seed = 333)

#mean of the (standard deviations of the  predictions per observation)
sd_cpuv <- apply(pred_cpuv, 2, sd) %>% mean
#sd roughly .07
#for the unconditional predictions from the multivariate model
sd_cp_mvun <- apply(pred_cp_gj_un$ychains[,1:140], 2, sd) %>% mean
#sd (.26) is a lot higher than for the univariate case
#for the conditional predictions
sd_cp_mvco <- apply(pred_cp_mvco$ychains[,1:140], 2, sd) %>% mean
#The sd (.25) is way higher than in the univariate case >> would mean that the conditional 
#predictions are a lot more uncertain; slightly lower than the unconditional predictions.
#So, the uncertainty is probably induced by the different modelling approaches (multi vs. 
#univariate + difference in fitting between gjam and rstanarm), not the conditioning. 

#Plotting uncertainties of predictions with boxplots for single randomly drawn observations
#draw the observation
i <- sample(seq_len(nrow(test)), size = 1)

#make the dataframe with all the predictions of observation i for all three prediction types
d_gg_bpso <- data.frame(c(pred_cpuv[,i], pred_cp_gj_un$ychains[,i], pred_cp_mvco$ychains[,i]))
names(d_gg_bpso) <- c("pred")
#adding a column with the prediction types
d_gg_bpso$type <- c(rep("Univariate", length(pred_cpuv[,i])),
                    rep("Unconditional Multivariate", length(pred_cp_gj_un$ychains[,i])),
                    rep("Conditional Multivariate", length(pred_cp_mvco$ychains[,i])))
#adding a column with the observation #
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

####For Anopheles####
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
#the sd (.22) is way higher than in the univariate case >> would mean that the conditional 
#predictions are a lot more uncertain; slightly lower than the unconditional predictions.
#Same discussion as for Culex perexiguus

#Plotting uncertainties of predictions with boxplots for single randomly drawn observations

#make the dataframe with all the predictions of observation i for all three prediction types
d_gg_bpso <- data.frame(c(pred_atuv[,i], pred_at_gj_un$ychains[,140+i], pred_at_mvco$ychains[,140+i]))
names(d_gg_bpso) <- c("pred")
#adding a column with the prediction types
d_gg_bpso$type <- c(rep("Univariate", length(pred_atuv[,i])),
                    rep("Unconditional Multivariate", length(pred_at_gj_un$ychains[,i])),
                    rep("Conditional Multivariate", length(pred_at_mvco$ychains[,i])))
#adding a column with observation #
d_gg_bpso$obs <- i
# Boxplot shows Boxes: 25% and 75 % Quantile; vertical line: median;
#lower whisker = smallest observation greater than or equal to lower hinge - 1.5 * IQR; 
#points: "remaining" outliers 
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
