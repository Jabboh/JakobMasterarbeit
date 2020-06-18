rm(list=ls())
setwd("C:\\Users\\jakob\\Documents\\JakobMasterarbeit")
#loading packages
library(readxl) #to read in the data
library(rstanarm) #Doing Bayesian probit GLMs for single species
library(pROC) #calculating the AUC
library(gjam) #Doing the joint estimation with the GJAM modelling aproach
library(ggplot2) #for plotting
library(dplyr) # For data manipulation
install.packages("matlab")
library(matlab)
####Data Preperation

#read in the data (Monthly species PA data for seven different mosquito species
#and according environmental covariates)
df <- read_excel("MonthlyData.xlsx")
str(df)
#We ignore An_atroparvus, because it does not seem to be PA-data.
#It looks like abundance data

#extract the PA data for all species
spec <- df[,7:14]
#deleting the An_atroparvus column
spec[,"An_atroparvus"] <- NULL
#To select four species for the analysis, we look for species that roughly occur in 50%
#of the observations.
summary(spec)
#Hence, we select Cxperpre, Cxpippre, Cxmodpre and Anatropre for the analysis:
y <- spec[,c("Cxpippre", "Cxmodpre", "Cxperpre", "Anatropre")]
#For these four species, Roiz et al. use the following covariates to explain their 
#monthly PAs: Inundation area (500 m buffer), NDVI (500 m buffer), Inundation area (2000 m),
#Inundation area before (200m) and NDVI month before (2000 m buffer). Since NDVI 500 m buffer
#and NDVI 2000 m buffer will be highly correlated, I will leave out NDVI 2000 m buffer.

#Normalizing Covariates
df[,17:ncol(df)] <- scale(df[,17:ncol(df)])

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

#fitting the model for Culex pipiens
fit_cpi <- stan_glm(Cxpippre ~ IA_500 + NDVI_500 + IABEF_2000 + NDVIBEF_2000, data = train, family = binomial(link = "probit"), seed = 333)
summary(fit_cpi)


#fitting the model for Culex modestus
fit_cm <- stan_glm(Cxmodpre ~ IA_500 + NDVI_500 + IABEF_2000 + NDVIBEF_2000, data = train, family = binomial(link = "probit"), seed = 333)
summary(fit_cm)

#fitting the model for Culex Perexiguus
fit_cp <- stan_glm(Cxperpre ~ IA_500 + NDVI_500 + IABEF_2000 + NDVIBEF_2000, data = train, family = binomial(link = "probit"), seed = 333)
summary(fit_cp)

#fitting the model for Anopheles atroparvus
fit_at <- stan_glm(Anatropre ~ IA_500 + NDVI_500 + IABEF_2000 + NDVIBEF_2000, data = train, family = binomial(link = "probit"), seed = 333)
summary(fit_at)

####predictions on the test set

#for Culex pipiens
pred_cpi_sin <- posterior_predict(fit_cpi, newdata = test, seed = 333)
pred_exp_cpi_sin <- colMeans((pred_cpi_sin))

#for Culex modestus
pred_cm_sin <- posterior_predict(fit_cm, newdata = test, seed = 333)
pred_exp_cm_sin <- colMeans((pred_cm_sin))

#for Culex perexiguus
pred_cp_sin <- posterior_predict(fit_cp, newdata = test, seed = 333)

#Is it correct that for each test data point 4000 draws from the posterior predictive 
#distribution are made? So I could take the average of these draws as an estimation of the 
# "expected" predicted y-value?
pred_exp_cp_sin <- colMeans((pred_cp_sin))

#for Anopheles atroparvus
pred_at_sin <- posterior_predict(fit_at, newdata = test, seed = 333)
pred_exp_at_sin <- colMeans((pred_at_sin))

####Quick Evaluation of Model Performance: To check whether model makes some sense

#For Culex pipiens
perf_cpi_sin <- auc(response = test$Cxpippre, predictor = pred_exp_cpi_sin)
perf_cpi_sin

#For Culex Modestus
perf_cm_sin <- auc(response = test$Cxmodpre, predictor = pred_exp_cm_sin)
perf_cm_sin

#For Culex perexiguus
perf_cp_sin <- auc(response = test$Cxperpre, predictor = pred_exp_cp_sin)
perf_cp_sin
#AUC = .74 >> seems ok for our purposes (remember: our goal is not to find the perfect
#model for our data, but rather evaluate whether knowing one species helps our predictions
#of the other one
perf_at_sin <- auc(response = test$Anatropre, predictor = pred_exp_at_sin)
perf_at_sin
#AUC = .79 >> seems OK for our purposes

####Fitting a joint model of both species with GJAM

#Define the model settings
types <- c("PA", "PA", "PA", "PA") #variable types of the responses
s <- length(types) # number of species
#y-data to train the model
y_train <- y[train_id, ]

#define model/algorithm parameters
ml   <- list(ng = 2000, burnin = 100, typeNames = types)
#is 2000 for number of Gibbs steps ok?

#runnig GJAM
joint <- gjam(~ IA_500 + NDVI_500 + IABEF_2000 + NDVIBEF_2000, ydata = y_train, xdata = train, modelList = ml)
summary(joint)

#Residual Correlations between the species
joint$parameters$corMu 
#is corMu actually the residual correlation? How do you calculate that?


####Comparison of the coefficients

#Make a table (dataframe) to store different coefficients and SEs according to the models
q = length(fit_cp$coefficients) # Number of predictors
#For Culex Perexiguus
#For the coefficients
cof_sum_px <- data.frame(matrix(ncol = 2, nrow = q))
colnames(cof_sum_px) <- c("Coefficients_sin", "Coefficients_gjam")
rownames(cof_sum_px) <- names(fit_cp$coefficients)

#For the SEs
se_sum_px <- data.frame(matrix(ncol = 2, nrow = q))
colnames(se_sum_px) <- c("SE_sin", "SE_gjam")
rownames(se_sum_px) <- names(fit_cp$coefficients)

###For the other species
#For the coefficients
cof_sum_at <- cof_sum_px
se_sum_at <- se_sum_px

cof_sum_cpi <- cof_sum_px
se_sum_cpi <- se_sum_px

cof_sum_cm <- cof_sum_px
se_sum_cm <- se_sum_px

#Filling the tables accordingly
cof_sum_px$Coefficients_sin <-  fit_cp$coefficients
cof_sum_px$Coefficients_gjam <- joint$parameters$betaMu[,"Cxperpre"]
cof_sum_at$Coefficients_sin <-  fit_at$coefficients
cof_sum_at$Coefficients_gjam <- joint$parameters$betaMu[,"Anatropre"]
cof_sum_cpi$Coefficients_sin <-  fit_cpi$coefficients
cof_sum_cpi$Coefficients_gjam <- joint$parameters$betaMu[,"Cxpippre"]
cof_sum_cm$Coefficients_sin <-  fit_cm$coefficients
cof_sum_cm$Coefficients_gjam <- joint$parameters$betaMu[,"Cxmodpre"]

se_sum_px$SE_sin <-  fit_cp$ses
se_sum_px$SE_gjam <- joint$parameters$betaSe[,"Cxperpre"]
se_sum_at$SE_sin <-  fit_at$ses
se_sum_at$SE_gjam <- joint$parameters$betaSe[,"Anatropre"]
se_sum_cpi$SE_sin <-  fit_cpi$ses
se_sum_cpi$SE_gjam <- joint$parameters$betaSe[,"Cxpippre"]
se_sum_cm$SE_sin <-  fit_cm$ses
se_sum_cm$SE_gjam <- joint$parameters$betaSe[,"Cxmodpre"]


#Coefficients and SEs for Culex perexiguus
cof_sum_px
se_sum_px
#Everything looks pretty similar, as we expected (We expexted the environmental coefficients to
#be the same. Difference between coefficients way smaller than according SEs. 

#For Culex pipiens
cof_sum_cpi
se_sum_cpi

#For Culex modestus
cof_sum_cm
se_sum_cm

#Coefficients and SEs for Anopheles troparvus
cof_sum_at
se_sum_at
#Everything looks pretty similar as expected. Difference between coefficients way smaller than according SEs.

####Conditional Prediction on the test set

# preparing testing data set for the gjam_predict function 
#y-data to test the model
y_test <- y[-train_id, ]
#x-data to test the model
test_gj <- df[-train_id,]
test_gj$intercept <- rep(1 , nrow(test_gj)) #adding the intercept
test_gj <- test_gj[,c("intercept", "IA_500", "NDVI_500", "IABEF_2000", "NDVIBEF_2000")] # getting rid of all the unused variables

#Doing the conditional prediction: Always one species conditioned on all the others
#Doing it for every species in a loop
for (i in names(y_test)){
  #else if-statement that gives us the appropriate abbreviation to work with in naming the variables
  if (i == "Cxpippre"){ 
    x = "cpi"
  } else if(i == "Cxmodpre"){
    x = "cm"
  } else if(i == "Cxperpre"){
    x = "cp"
  } else {
    x = "at"
  }
  #define newdata so that one species is left out (the one that will be predicted)
  newdata <- list(xdata = test_gj, ydataCond = select(y_test, -i), nsim = 200) # conditionally predict out-of-sample
  assign(paste0("p_", x), gjamPredict(output = joint, newdata = newdata))
}


###Model Performance Evaluation with AUC

#For Culex pipiens
perf_cpi_gj <- auc(response = test$Cxpippre, predictor = p_cpi$prPresent[,"Cxpippre"])
perf_cpi_gj
# The AUC is with .81 considerably better than for the univariate case (.63)

#For Culex modestus
perf_cm_gj <- auc(response = test$Cxmodpre, predictor = p_cm$prPresent[,"Cxmodpre"])
perf_cm_gj
# The AUC is with .68 slightly better than for the univariate case (.67)

#For Culex Perexiguus
perf_cp_gj <- auc(response = test$Cxperpre, predictor = p_cp$prPresent[,"Cxperpre"])
perf_cp_gj
# The AUC is with .83 considerably better than for the univariate case (.74)

#For Anopheles atroparvus
perf_at_gj <- auc(response = test$Anatropre, predictor = p_at$prPresent[,"Anatropre"])
perf_at_gj
#The AUC is with .82 slightly better than for the univariate case (.79)

#Overall, the improvement is stronger than for the case of only 2 species. But, there is considerable heterogeneity.
#It seems like for species with low residual correlations with other species (here: Culex modestus) and species for
#which environmental variables already explain a lot (here: Anopheles atroparvus), the improvement in AUC
#is very low.

####plot GJAM-Predictions against the univariate probit predictions

#for Culex perexiguus
#make a categorical variable for which species are present: 0: no species present; 1: cpi, 2: cm; 3:at
#4:cpi & cm, 5: cpi & at; 6: cm & at; 7: all are present. This will later be used for coloring in plot
#select all the species on which we condition
pres_cp <- select(as.data.frame(p_cp$prPresent), -"Cxperpre")
#Function that returns categorical value corresponding to the above specification (value shows which species
#are present)
presence <- function (x){
  if(x[1] == 0 & x[2] == 0 & x[3] == 0){
    zeros = 0
  } else if(x[1] == 1 & x[2] == 0 & x[3] == 0){
    zeros = 1
  } else if(x[1] == 0 & x[2] == 1 & x[3] == 0){
    zeros = 2
  } else if(x[1] == 0 & x[2] == 0 & x[3] == 1){
    zeros = 3
  } else if(x[1] == 1 & x[2] == 1 & x[3] == 0){
    zeros = 4
  } else if(x[1] == 1 & x[2] == 0 & x[3] == 1){
    zeros = 5
  } else if(x[1] == 0 & x[2] == 1 & x[3] == 1){
    zeros = 6
  } else {
    zeros = 7 
  } 
  return(zeros)
}
#apply function to all rows of dataframe that shows PA of all species we condition on (So, basically
#we transofrm a three dimensional variable into a one-dimensional one without losing any information.
#This is needed for ggplot later)
pres_cat_cp <- apply(pres_cp, 1, presence)


d_gg_cp <- data.frame(cbind(pred_exp_cp_sin, p_cp$prPresent[,"Cxperpre"],  pres_cat_cp, y_test$Cxperpre))
names(d_gg_cp) <- c("cp_pr", "cp_gj", "con_spec", "cp")
provsgj_cp <- ggplot(d_gg_cp, aes(x=cp_pr, y=cp_gj) ) +
  geom_point(aes(fill = factor(con_spec), shape = factor(cp)), color= "black") +
  scale_shape_manual(values=c(21, 22), labels = c("Absent", "Present")) + 
  scale_fill_manual(values=c("white", "red", "blue", "yellow", "purple", "orange", "green", "black"),
                    labels = c("all absent", "cpi", "cm", "at", "cpi & cm", "cpi & at", "cm & at", "all present")) +
  ggtitle("Univariate vs. Conditional Predictions for Culex Perexiguus") + 
  xlab("Predictions from Univariate Probit ") + ylab("Conditional Predictions from GJAM") +
  guides(fill=guide_legend(override.aes=list(shape=21))) +
  labs(fill = "Presence of Conditioning Species", shape = "PA of Culex Perexiguus") #Conditioning Species kann man bestimmt nicht sagen

provsgj_cp

#You can see that there is a positive linear relationship between the predictions
#of the two models. This indicates that both models roughly do the same/environmental signals are 
#treated similarily.You can see that there are heterogeneous effects
#of the presences of the different species we condition on. The presence of Culex Pipiens and Anopheles
#troparvus has a positive effect on the predictions in GJAM, whereas the mere presence of Culex modestus has
#a negative effect. There seem to be three lines around which the points gather: 1. The bottom line:
#Absence of all species and presence of just Culex modestus. (Probit predictions are higher than GJAM
#predictions); 2. The middle line (GJAM prediction are slightly higher than probit predictions): Presence of (i) Culex pipiens with or without Culex Modestus , (ii)
#Anopheles Atroparvus with or without Culex modestus; 3. Upper line(GJAM predictions are much higher than
#probit predictions): Culex pipiens and Anopheles troparvus are present with or without Culex modestus.
#This, can be explained with the fact that residual correlations of Culex Perexiguus are high with 
#Culex pipiens and Anopheles Atroparvus, but very low with Culex modestus (the raw correlation between
#species has the same tendencies)


###for Anopheles Atroparvus
#Dataframe with species that we conditioned on:
pres_at <- select(as.data.frame(p_at$prPresent), -"Anatropre")

#Make one categorical variable out of above datafame with the following levels:
#0: no species present; 1: cpi, 2: cm; 3:cp 4:cpi & cm, 5: cpi & cp; 6: cm & cp; 7: all are present. 
pres_cat_at <- apply(pres_at, 1, presence)

#Input data for ggplot
d_gg_at <- data.frame(cbind(pred_exp_at_sin, p_at$prPresent[,"Anatropre"],  pres_cat_at, y_test$Anatropre))
names(d_gg_at) <- c("at_pr", "at_gj", "con_spec", "at")
provsgj_at <- ggplot(d_gg_at, aes(x=at_pr, y=at_gj) ) +
  geom_point(aes(fill = factor(con_spec), shape = factor(at)), color= "black") +
  scale_shape_manual(values=c(21, 22), labels = c("Absent", "Present")) + 
  scale_fill_manual(values=c("white", "red", "blue", "yellow", "purple", "orange", "green", "black"),
                    labels = c("all absent", "cpi", "cm", "cp", "cpi & cm", "cpi & cp", "cm & cp", "all present")) +
  ggtitle("Univariate vs. Conditional Predictions for Anopheles Atroparvus") + 
  xlab("Predictions from Univariate Probit ") + ylab("Conditional Predictions from GJAM") +
  guides(fill=guide_legend(override.aes=list(shape=21))) +
  labs(fill = "Presence of Conditioning Species", shape = "PA of Anopheles Atroparvus")
provsgj_at
#Similar results as for Culex Perexiguus.

###For Culex Pipiens

#Dataframe with species that we conditioned on:
pres_cpi <- select(as.data.frame(p_cpi$prPresent), -"Cxpippre")

#Make one categorical variable out of above datafame with the following levels:
#0: no species present; 1: cm, 2: cp; 3:at 4:cm & cp, 5: cm & at; 6: cp & at; 7: all are present. 
pres_cat_cpi <- apply(pres_cpi, 1, presence)

#Input data for ggplot
d_gg_cpi <- data.frame(cbind(pred_exp_cpi_sin, p_cpi$prPresent[,"Cxpippre"],  pres_cat_cpi, y_test$Cxpippre))
names(d_gg_cpi) <- c("cpi_pr", "cpi_gj", "con_spec", "cpi")
provsgj_cpi <- ggplot(d_gg_cpi, aes(x=cpi_pr, y=cpi_gj) ) +
  geom_point(aes(fill = factor(con_spec), shape = factor(cpi)), color= "black") +
  scale_shape_manual(values=c(21, 22), labels = c("Absent", "Present")) + 
  scale_fill_manual(values=c("white", "red", "blue", "yellow", "purple", "orange", "green", "black"),
                    labels = c("all absent", "cm", "cp", "at", "cm & cp", "cm & at", "cp & at", "all present")) +
  ggtitle("Univariate vs. Conditional Predictions for Culex Pipiens") + 
  xlab("Predictions from Univariate Probit ") + ylab("Conditional Predictions from GJAM") +
  guides(fill=guide_legend(override.aes=list(shape=21))) +
  labs(fill = "Presence of Conditioning Species", shape = "PA of Culex Modestus")
provsgj_cpi

#Linear relationships are hardly detectable (Why is this the case?). But you can see which species lead to higher predictions in
#gjam: Presence of Culex perexiguus and Anopheles Atroparvus lead to highest increase; followed by 
#Culex Modestus and Culex perexiguus as well as Culex Modestus and Anopheles Atroparvus; followed
#by Culex perexiguus as well as Anopheles Atroparvus. GJAM predictions are only consistently higher 
#than univariate probit predictions if Culex perexiguus and Anopheles Atroparvus are both present.
#Residual and raw correlation, again, explain which species leads to higher predictions in GJAM.
###For Culex modestus

#Dataframe with species that we conditioned on:
pres_cm <- select(as.data.frame(p_cm$prPresent), -"Cxmodpre")

#Make one categorical variable out of above datafame with the following levels:
#0: no species present; 1: cpi, 2: cp; 3:at 4:cpi & cp, 5: cpi & at; 6: cp & at; 7: all are present. 
pres_cat_cm <- apply(pres_cm, 1, presence)

#Input data for ggplot
d_gg_cm <- data.frame(cbind(pred_exp_cm_sin, p_cm$prPresent[,"Cxmodpre"],  pres_cat_cm, y_test$Cxmodpre))
names(d_gg_cm) <- c("cm_pr", "cm_gj", "con_spec", "cm")
provsgj_cm <- ggplot(d_gg_cm, aes(x=cm_pr, y=cm_gj) ) +
  geom_point(aes(fill = factor(con_spec), shape = factor(cm)), color= "black") +
  scale_shape_manual(values=c(21, 22), labels = c("Absent", "Present")) + 
  scale_fill_manual(values=c("white", "red", "blue", "yellow", "purple", "orange", "green", "black"),
                    labels = c("all absent", "cpi", "cp", "at", "cpi & cp", "cpi & at", "cp & at", "all present")) +
  ggtitle("Univariate vs. Conditional Predictions for Culex Modestus") + 
  xlab("Predictions from Univariate Probit ") + ylab("Conditional Predictions from GJAM") +
  guides(fill=guide_legend(override.aes=list(shape=21))) +
  labs(fill = "Presence of Conditioning Species", shape = "PA of Culex Modestus")
provsgj_cm
#Pretty much looks like one "Punktewolke" (because residual and raw correlations are very weak),
#only very small difference between the "conditioning" species. 
#The predictions between GJAM and Probit seem to be pretty similar across observations (consistent with correlations).
#GJAM predicts slightly higher presence probabilities, if Culex Pipiens and or Anopheles Atroparvus
#are present and slightly lower probabilities, if all species are absent or only Culex Perexigus is present.

#######Unconditional Predictions

###Unconditional GJAM Predictions 
newdata <- list(xdata = test_gj, nsim = 200) # conditionally predict out-of-sample
#Doing the actual prediction
p_uc      <- gjamPredict(output = joint, newdata = newdata)

###AUCs
#For Culex Pipiens
perf_cpi_uc <- auc(response = test$Cxpippre, predictor = p_uc$prPresent[,"Cxpippre"])
perf_cpi_uc
# The AUC is with .64 pretty much the same as in the univariate case (.63)

# For Culex Modestus
perf_cm_uc <- auc(response = test$Cxmodpre, predictor = p_uc$prPresent[,"Cxmodpre"])
perf_cm_uc
# The AUC is with .66 pretty much the same as in the univariate case (.67)

#For Culex Perexiguus
perf_cp_uc <- auc(response = test$Cxperpre, predictor = p_uc$prPresent[,"Cxperpre"])
perf_cp_uc
# The AUC is with .74 pretty much the same as in the univariate case (.74)

#For Anopheles atroparvus
perf_at_uc <- auc(response = test$Anatropre, predictor = p_uc$prPresent[,"Anatropre"])
perf_at_uc
#The AUC is with .79 pretty much the same as in the univariate case (.79)

#Overall, there are no improvements in GJAM compared to a univariate probit model.

### plot unconditional GJAM-predictions against the univariate prediction

#For Culex Perexiguus
d_gg_uc_cp <- data.frame(cbind(pred_exp_cp_sin, p_uc$prPresent[,"Cxperpre"], y_test$Cxperpre))
names(d_gg_uc_cp) <- c("cp_pr", "cp_gj_un", "cp")
provsgjun_cp <- ggplot(d_gg_uc_cp, aes(x=cp_pr, y=cp_gj_un, color=factor(cp))) + geom_point() +
  ggtitle("Univariate vs. Unconditional GJAM Predictions for Culex Perexiguus") +
  xlab("Predictions from Univariate Probit ") + ylab("Unconditional Predictions from GJAM") +
  labs(color = "True PA of Culex Perexiguus")
provsgjun_cp
#They look as if they are centered around the identity line >> the predictions are more or less the same!

#for Anopheles Atroparvus
d_gg_uc_at <- data.frame(cbind(pred_exp_at_sin, p_uc$prPresent[,"Anatropre"], y_test$Anatropre))
names(d_gg_uc_at) <- c("at_pr", "at_gj_un", "at")
provsgjun_at <- ggplot(d_gg_uc_at, aes(x=at_pr, y=at_gj_un, color=factor(at))) + geom_point() +
  ggtitle("Univariate vs. Unconditional GJAM Predictions for Anopheles Atroparvus") +
  xlab("Predictions from Univariate Probit ") + ylab("Unconditional Predictions from GJAM") +
  labs(color = "True PA of Anopheles Atroparvus")
provsgjun_at
#They look as if they are centered around the identity line >> the predictions are more or less the same!
#I think that prooves Björn Point >> We could also make one plot with all the species. Since they should be
#all on the same line, this could illustrate the finding, that unconditional GJAM predictions and univariate
#predictions are pretty much the same.
#In summary, if we would plot conditional GJAM predictions against unconditional gjam predictions,
#we would get the same results as in our first two plots (probit predictions vs. conditional gjam predictions)

