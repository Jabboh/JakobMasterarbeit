#stan_glm function macht Faxen, der "initialization  failed"-error macht die Probleme

rm(list=ls())
setwd("C:\\Users\\jakob\\Documents\\JakobMasterarbeit\\Data")
#loading packages
library(readxl) #to read in the data
library(rstanarm) #Doing Bayesian probit GLMs for single species
library(pROC) #calculating the AUC
library(gjam) #Doing the joint estimation with the GJAM modelling aproach
library(ggplot2) #for plotting
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
#To select two species for a first analysis, we look for species that roughly occur in 50%
#of the observations.
summary(spec)
#Hence, we select Cxperpre and Anatropre for our first analysis:
y <- spec[,c("Cxperpre", "Anatropre")]
#For these two species, Roiz et al. use the following covariates to explain their 
#monthly PAs: Inundation area (500 m buffer), NDVI (500 m buffer), Inundation area (2000 m)
#and NDVI month before (2000 m buffer). Since NDVI 500 m buffer and NDVI 2000 m buffer
#will be highly correlated, I will leave out NDVI 2000 m buffer.

#Normalizing Covariates: Otherwise, interaction terms hardy interpretable and skewed
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

#fitting the model for Culex Perexiguus
fit_cp <- stan_glm(Cxperpre ~ (IA_500 + NDVI_500 + NDVIBEF_2000)^2 + I(IA_500^2) + I(NDVI_500^2) + I(NDVIBEF_2000^2), data = train, family = binomial(link = "probit"), init_r = 1, seed = 333)
summary(fit_cp)
#"significant on the .1-level are IA_500, NDVIBEF_2000, I(NDVI_500^2), I(NDVIBEF_2000^2), IA_500:NDVIBEF_2000
#and NDVI_500:NDVIBEF_2000
#Remember without the interactions and quadratic terms the following variables were significant:
#(Intercept), IA_500 and NDVIBEF_2000

#fitting the model for Anopheles atroparvus
fit_at <- stan_glm(Anatropre ~ (IA_500 + NDVI_500 + NDVIBEF_2000)^2 + I(IA_500^2) + I(NDVI_500^2) + I(NDVIBEF_2000^2), data = train, family = binomial(link = "probit"), init_r = 1, seed = 333)
summary(fit_at)
#"significant on the .1-level are NDVIBEF_2000, I(NDVIBEF_2000^2), IA_500:NDVI_500 and IA_500:NDVIBEF_2000
#Remember without the interactions and quadratic terms the following variables were significant: Intercept
#and NDVIBEF_2000
####predictions on the test set

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
perf_cp_sin <- auc(response = test$Cxperpre, predictor = pred_exp_cp_sin)
perf_cp_sin
#AUC = .73. This is the same as in the case without the quadratic and interaction terms >> Does that make sense?
#should be   a little higher, bc way more predictors (but don't forget overfitting!...)

perf_at_sin <- auc(response = test$Anatropre, predictor = pred_exp_at_sin)
perf_at_sin
#AUC = .8 (.03 higher than in our "simple" model)

####Fitting a joint model of both species with GJAM

#Define the model settings
types <- c("PA", "PA") #variable types of the responses
s <- length(types) # number of species
#y-data to train the model
y_train <- y[train_id, ]

#define model/algorithm parameters
ml   <- list(ng = 2000, burnin = 100, typeNames = types)
#is 2000 for number of Gibbs steps ok?

#runnig GJAM
joint <- gjam(~ (IA_500 + NDVI_500 + NDVIBEF_2000)^2 + I(IA_500^2) + I(NDVI_500^2) + I(NDVIBEF_2000^2), ydata = y_train, xdata = train, modelList = ml)
summary(joint)

#Residual Correlations between the species
joint$parameters$corMu 
#is corMu actually the residual correlation? How do you calculate that?
#You first estimate the environmental parameters and then you produce fitted y-values with the coefficients.
#2. You take the residuals of this estimation and calculate the correlation between the residuals
#of the two speices?
cor(train$Cxperpre, train$Anatropre)
#Isn't it improbable that the residual correlation is bigger than the "raw" correlation between the species, assuming
#that both species favor the same environmental niches?

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

###For Anopheles troparvus
#For the coefficients
cof_sum_at <- cof_sum_px
se_sum_at <- se_sum_px

#Filling the tables accordingly
cof_sum_px$Coefficients_sin <-  fit_cp$coefficients
cof_sum_px$Coefficients_gjam <- joint$parameters$betaMu[,"Cxperpre"]
cof_sum_at$Coefficients_sin <-  fit_at$coefficients
cof_sum_at$Coefficients_gjam <- joint$parameters$betaMu[,"Anatropre"]

se_sum_px$SE_sin <-  fit_cp$ses
se_sum_px$SE_gjam <- joint$parameters$betaSe[,"Cxperpre"]
se_sum_at$SE_sin <-  fit_at$ses
se_sum_at$SE_gjam <- joint$parameters$betaSe[,"Anatropre"]

#Coefficients and SEs for Culex perexiguus
cof_sum_px
se_sum_px
#Everything looks pretty similar, as we expected (We expexted the environmental coefficients to
#be the same. Difference between coefficients way smaller than according SEs. 

#Coefficients and SEs for Anopheles troparvus
cof_sum_at
se_sum_at
#Everything looks pretty similar as expected. Difference between coefficients way smaller than according SEs.

####Conditional Prediction on the test set

# preparing testing data set for the gjam_predict function 

#adding all the squared, interaction and intercept terms (This is cumbersome, but I don't know a way around
#in "gjamPredict")

df$"I(IA_500^2)" <- (df$IA_500)^2
df$"I(NDVI_500^2)" <- (df$NDVI_500)^2
df$"I(NDVIBEF_2000^2)" <- (df$NDVIBEF_2000)^2
df$"IA_500:NDVI_500" <- (df$IA_500) * (df$NDVI_500)
df$"IA_500:NDVIBEF_2000" <- (df$IA_500) * (df$NDVIBEF_2000)
df$"NDVI_500:NDVIBEF_2000" <- (df$NDVIBEF_2000) * (df$NDVI_500)


#rerun gjam with this structure
joint2 <- gjam(~ IA_500 + NDVI_500 + NDVIBEF_2000 + I(IA_500^2) + I(NDVI_500^2) + I(NDVIBEF_2000^2) +
                 IA_500:NDVI_500 + IA_500:NDVIBEF_2000 + NDVI_500:NDVIBEF_2000, ydata = y_train,
               xdata = train, modelList = ml)


#y-data to test the model
y_test <- y[-train_id, ]
#x-data to test the model
test_gj <- df[-train_id,]
test_gj$intercept <- rep(1 , nrow(test_gj)) #adding the intercept

test_gj <- test_gj[,c("intercept", "IA_500", "NDVI_500", "NDVIBEF_2000", "I(IA_500^2)",
                      "I(NDVI_500^2)", "I(NDVIBEF_2000^2)", "IA_500:NDVI_500", "IA_500:NDVIBEF_2000",
                      "NDVI_500:NDVIBEF_2000")] # getting rid of all the unused variables

#Culex perexiguus conditioned on Anopheles atroparvus
#storing input data in newdata
newdata <- list(xdata = test_gj, ydataCond = y_test[,2], nsim = 200) # conditionally predict out-of-sample
#Doing the actual prediction
p_cp      <- gjamPredict(output = joint2, newdata = newdata)

#Anopheles atroparvus conditioned on Culex perexiguus 
#preparing newdata
newdata <- list(xdata = test_gj, ydataCond = y_test[,1], nsim = 200) # conditionally predict out-of-sample
#Doing the actual prediction
p_at     <- gjamPredict(output = joint2, newdata = newdata)
###Model Performance Evaluation with AUC
#For Culex Perexiguus
perf_cp_gj <- auc(response = test$Cxperpre, predictor = p_cp$prPresent[,1])
perf_cp_gj
# The AUC is with .77 slightly worse than for the simple model (.78)

#For Anopheles atroparvus
perf_at_gj <- auc(response = test$Anatropre, predictor = p_at$prPresent[,2])
perf_at_gj
#The AUC is with .8 slightly better than for the simple model (.79)
#Overall, there are improvements.

####plot GJAM-Predictions against the univariate probit predictions

#for Culex perexiguus
d_gg_cp <- data.frame(cbind(pred_exp_cp_sin, p_cp$prPresent[,1], p_cp$prPresent[,2], y_test$Cxperpre))
names(d_gg_cp) <- c("cp_pr", "cp_gj", "at", "cp")
provsgj_cp <- ggplot(d_gg_cp, aes(x=cp_pr, y=cp_gj, color=factor(at), shape = factor(cp))) + geom_point() +
  ggtitle("Univariate vs. Conditional Predictions for Culex Perexiguus") +
  xlab("Predictions from Univariate Probit ") + ylab("Conditional Predictions from GJAM") +
  labs( color = "PA of Anopheles Atroparvus", shape = "PA of Culex Perexiguus")
provsgj_cp
#It looks roughly the same as in the simple model, but details are different of distribution are different

#for Anopheles Atroparvus
d_gg_at <- data.frame(cbind(pred_exp_at_sin, p_at$prPresent[,2], p_at$prPresent[,1], y_test$Anatropre))
names(d_gg_at) <- c("at_pr", "at_gj", "cp", "at")
provsgj_at <- ggplot(d_gg_at, aes(x = at_pr, y = at_gj, color=factor(cp), shape = factor(at))) + geom_point() +
  ggtitle("Univariate vs. Conditional Predictions for Anopheles Atroparvus") +
  xlab("Predictions from Univariate Probit ") + ylab("Conditional Predictions from GJAM") +
  labs( color = "PA of Culex Perexiguus", shape = "PA of Anopheles Atroparvus")
provsgj_at
#Similar results as for Culex Perexiguus.

