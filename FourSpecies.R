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
#4:cpi & cm, 5: cm & at; 6: at & cpi; 7: all are present
pres <- select(as.data.frame(p_cp$prPresent), -"Cxperpre")
function (x){
  
}
pres_cat <- zeros(1, 1)
for (i in ncol(pres)){
  if(pres[1,i] == 0 & pres[2,i] == 0 & pres[3,i] == 0){
    zeros = 0
  } else if(pres[1,i] == 1 & pres[2,i] == 0 & pres[3,i] == 0){
    zeros = 1
  } else if(pres[1,i] == 0 & pres[2,i] == 1 & pres[3,i] == 0){
    zeros = 2
  } else if(pres[1,i] == 0 & pres[2,i] == 0 & pres[3,i] == 1){
    zeros = 3
  } else if(pres[1,i] == 1 & pres[2,i] == 1 & pres[3,i] == 0){
    zeros = 4
  } else if(pres[1,i] == 1 & pres[2,i] == 0 & pres[3,i] == 1){
    zeros = 5
  } else if(pres[1,i] == 0 & pres[2,i] == 1 & pres[3,i] == 1){
    zeros = 6
  } else if(pres[1,i] == 1 & pres[2,i] == 1 & pres[3,i] == 1){
    zeros = 7
  }
}
pres[1,1]
d_gg_cp <- data.frame(cbind(pred_exp_cp_sin, p_cp$prPresent[,"Cxperpre"],  select(as.data.frame(p_cp$prPresent), -"Cxperpre"), y_test$Cxperpre))
names(d_gg_cp) <- c("cp_pr", "cp_gj", "cpi", "cm" ,"at", "cp")
provsgj_cp <- ggplot(d_gg_cp, aes(x=cp_pr, y=cp_gj, color=factor("cpi"), fill = factor("cm"), shape = factor(cp))) +
  geom_point() + ggtitle("Univariate vs. Conditional Predictions for Culex Perexiguus") +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  xlab("Predictions from Univariate Probit ") + ylab("Conditional Predictions from GJAM") +
  labs( color = "PA of Anopheles Atroparvus", shape = "PA of Culex Perexiguus")
provsgj_cp
#You can see that there is a positive linear relationship between the predictions
#of the two models grouped by the PA of Anopheles Atroparvus (the species we conditioned on). This indicates
#that both models roughly do the same/environmental signals are treated similarily.
#The slopes of the two lines are not equal to 1 (Would we expect this? I think, we do, if we 
#assume that the environmental coefficients are the same). You can see that the conditioning on 
#Anapheles Atroparvus has a clear effect (The blue and red points form two distinct groups)
#on the predictions in GJAM compared to the univariate predictions. GJAM predicts a roughly 37 % points
#higher probability of occurence for plots where Anopheles Atroparvus is present compared to
#plots where it's absent.

#for Anopheles Atroparvus
d_gg_at <- data.frame(cbind(pred_exp_at_sin, p_at$prPresent[,2], p_at$prPresent[,1], y_test$Anatropre))
names(d_gg_at) <- c("at_pr", "at_gj", "cp", "at")
provsgj_at <- ggplot(d_gg_at, aes(x = at_pr, y = at_gj, color=factor(cp), shape = factor(at))) + geom_point() +
  ggtitle("Univariate vs. Conditional Predictions for Anopheles Atroparvus") +
  xlab("Predictions from Univariate Probit ") + ylab("Conditional Predictions from GJAM") +
  labs( color = "PA of Culex Perexiguus", shape = "PA of Anopheles Atroparvus")
provsgj_at
#Similar results as for Culex Perexiguus.



grp = rep(LETTERS[1:5],each = 2)