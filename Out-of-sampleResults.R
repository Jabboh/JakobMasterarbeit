#############################Results of my analyis 2: Out-of Sample
#1. Data preparation and fitting of final model
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
library(corrplot)
library(loo) #to calculate WAIC
library(bayesplot) #Some handy features for plotting in the Bayesian realm
library(matrixStats)#for quantile calculations
library(ggplot2) #for plotting

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
##################################################6.2. Out-of-Sample
#a. Conditional Predictions of gjam vs. Unconditional Predictions of 
#gjam vs. predictions of univariate model

#We evaluate predictive performance on our test set with the AUC
#####For Culex perexiguus
#We make the three predictions ((i) univariate, (ii) unconditional multivariate and (iii)
#conditional multivariate)

###(i) univariate predictions on the observation scale
pred_cp_uv <- posterior_predict(fit_fin_cp, newdata = test, seed = 333)
#take the average of these draws per observation as an estimation of the 
# "expected" predicted y-value/predicted probability of occurence
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
#make a factor out cp. This comes in handy for the plotting
d_gg_uc_cp$cp <- factor(d_gg_uc_cp$cp, levels = c(0, 1), labels = c("absent", "present"))

uvvsmvun_cp <- ggplot(d_gg_uc_cp, aes(x=cp_uv, y=cp_mv_un, shape=factor(cp))) + geom_point() +
  ggtitle("Univariate vs. Unconditional Multivariate Predictions for Culex perexiguus") +
  xlab("Predictions from Univariate Probit ") + 
  ylab("Unconditional Predictions from Multivariate Probit") +
  labs(shape = "Observed PA of Culex perexiguus") + ylim(0, 1)
uvvsmvun_cp
#They look as if they are centered around the identity line >> the predictions are more or less 
#the same! There is a slight tendency, though: For low predictions multivariate predictions are
#higher than univariate, but for high predictions univariate predictions are higher than 
#multivariate ones. This weakly contradicts hypothesis 3?

###Conditional predictions
#Culex perexiguus conditioned on Anopheles troparvus
#storing input data in newdata, inclduing the PA of Anopheles
set.seed(333)
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

#Plotting Univariate Predictions vs. Conditional Multivariate Predictions
d_gg_cp <- data.frame(cbind(pred_cp_uv, p_cp_mvco, pred_cp_mvco$prPresent[,2], y_test$Cxperpre))
names(d_gg_cp) <- c("cp_uv", "cp_mv", "at", "cp")
#Make the PA-variables to factors
d_gg_cp$cp <- factor(d_gg_cp$cp, levels = c(0, 1), labels = c("absent", "present"))
d_gg_cp$at <- factor(d_gg_cp$at, levels = c(0, 1), labels = c("absent", "present"))

provsgj_cp <- ggplot(d_gg_cp, aes(x=cp_uv, y=cp_mv, color=at, shape = cp)) + geom_point() +
  ggtitle("Univariate vs. Conditional Multivariate Predictions for Culex perexiguus") +
  xlab("Predictions from Univariate Probit ") + ylab("Conditional Predictions from Multivariate Probit") +
  labs( color = "PA of Anopheles troparvus", shape = "PA of Culex perexiguus") + ylim(0,1)
provsgj_cp
#You can see that there is a positive linear relationship between the predictions
#of the two models grouped by the PA of Anopheles troparvus (the species we conditioned on).
#This indicates
#that both models roughly do the same/environmental signals are treated similarily.
#The slopes of the two lines are not equal to 1 (Would we expect this? I think, we do, if we 
#assume that the environmental coefficients are the same PUH, I still dont know the answer
#to this question). You can see that the conditioning on 
#Anapheles troparvus has a clear effect (The blue and red points form two distinct groups)
#on the multivariate predictions compared to the univariate predictions. The multivariate model
#predicts a roughly 12.5 % points
#higher probability of occurence for plots where Anopheles troparvus is present compared to
#plots where it's absent.



#####for Anopheles troparvus
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
d_gg_uc_at <- data.frame(cbind(pred_at_uv, p_at_gj_un, y_test$Anatropre))
names(d_gg_uc_at) <- c("at_uv", "at_mv_un", "at")
d_gg_uc_at$at <- factor(d_gg_uc_at$at, levels = c(0, 1), labels = c("absent", "present"))

uvvsmvun_at <- ggplot(d_gg_uc_at, aes(x=at_uv, y=at_mv_un, shape=factor(at))) + geom_point() +
  ggtitle("Univariate vs. Unconditional Multivariate Predictions for Culex perexiguus") +
  xlab("Predictions from Univariate Probit ") + 
  ylab("Unconditional Predictions from Multivariate Probit") +
  labs(shape = "True PA of Culex perexiguus") + ylim(0, 1)
uvvsmvun_at
#Predictions are on one line, but not the identity line >> contradicts our hypothesis that 
#univariate probit predictions and unconditional multivariate probit predictions do not differ
#As for Culex, mv predictions are higher for low predictions and the other way around for high
#predictions

###Conditional predictions
#Anopheles troparvus conditioned on Culex perexiguus 
#storing input data in newdata, inclduing the PA of Culex
set.seed(333)
newdata <- list(xdata = test_gj, ydataCond = y_test[,1], nsim = 4000) # conditionally predict out-of-sample
#Doing the actual prediction
pred_at_mvco      <- gjamPredict(output = joint_fin, newdata = newdata)
#retrieving the predictions (the mean of the simulations)
p_at_mvco <- pred_at_mvco$sdList$yMu[,2] 
#AUC
auc_at_mvco <- auc(response = test$Anatropre, predictor = p_at_mvco)
auc_at_mvco
#AUC roughly .89 >> So, there is a considerable improvement of .04 in AUC >> corroborates
#our hypothesis that PA data enhances predictions :)

#Plotting Univariate Predictions vs. Conditional Multivariate Predictions
d_gg_at <- data.frame(cbind(pred_at_uv, p_at_mvco, y_test$Anatropre, y_test$Cxperpre))
names(d_gg_at) <- c("at_uv", "at_mv", "at", "cp")
#Make the PA-variables to factors
d_gg_at$cp <- factor(d_gg_at$cp, levels = c(0, 1), labels = c("absent", "present"))
d_gg_at$at <- factor(d_gg_at$at, levels = c(0, 1), labels = c("absent", "present"))

provsgj_at <- ggplot(d_gg_at, aes(x=at_uv, y=at_mv, color=cp, shape = at)) + geom_point() +
  ggtitle("Univariate vs. Conditional Multivariate Predictions for Anopheles troparvus") +
  xlab("Predictions from Univariate Probit ") + ylab("Conditional Predictions from Multivariate Probit") +
  labs( color = "PA of Culex perexiguus", shape = "PA of Anopheles troparvus") + ylim(0,1)
provsgj_at

#We can see that the conditioning on the PA of CUlex has a clear effect on predictions of 
#mv-probit. It's roughly .125 high, if Culex is present. Moreover, we detect the same pattern
#as in the unconditional vs. uv-plot: for low predictions mv has higher predictions and for 
#high predictions the other way around.

####################################6.2b Uncertainty in the predictions
#We take a look at the dispersion of the predictions for the three prediction types.  

######For Culex perexiguus
###For the univariate model
pred_cpuv <- posterior_linpred(fit_fin_cp, transform = T, newdata = test, seed = 333)

#mean of the (standard deviations of the  predictions per observation)
sd_cpuv <- apply(pred_cpuv, 2, sd) %>% mean
#sd roughly .07
#for the unconditional predictions
sd_cp_mvun <- apply(pred_cp_gj_un$ychains[,1:140], 2, sd) %>% mean
#sd (.26) is a lot higher than for the univariate case
#for the conditional predictions
sd_cp_mvco <- apply(pred_cp_mvco$ychains[,1:140], 2, sd) %>% mean
#DAMN, the sd (.25) is way higher than in the univariate case >> would mean that the condtional 
#predictions are a lot more uncertain; slightly lower than the unconditional predictions.
#So, the uncertainty is probably induced by the different modelling approaches (multi vs. 
#univariate + difference in fitting between gjam and rstanarm), not the conditioning. I think,
#the reason is that gjam additonally (to the linear predictor with the environmental 
#covariates)draws the residual error from the multivariate normal distribution>> 
#More dispersion

#Plotting uncertainties of predictions with boxplots for single randomly drawn observations
i <- sample(seq_len(nrow(test)), size = 1)

#make the dataframe with all the predictions of observation i for all three prediction types
d_gg_bpso <- data.frame(c(pred_cpuv[,i], pred_cp_gj_un$ychains[,i], pred_cp_mvco$ychains[,i]))
names(d_gg_bpso) <- c("pred")
#adding a column with the prediction types
d_gg_bpso$type <- c(rep("Univariate", length(pred_cpuv[,i])),
                    rep("Unconditional Multivariate", length(pred_cp_gj_un$ychains[,i])),
                    rep("Conditional Multivariate", length(pred_cp_mvco$ychains[,i])))
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

####################For Anopheles
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
#the sd (.22) is way higher than in the univariate case >> would mean that the condtional 
#predictions are a lot more uncertain; slightly lower than the unconditional predictions.
#So, the uncertainty is probably induced by the different modelling approaches (multi vs. 
#univariate + difference in fitting between gjam and rstanarm), not the conditioning. I think,
#the reason is that gjam additonally (to the linear predictor with the environmental 
#covariates)draws the residual error from the multivariate normal distribution>> 
#More dispersion

#Plotting uncertainties of predictions with boxplots for single randomly drawn observations

#make the dataframe with all the predictions of observation i for all three prediction types
d_gg_bpso <- data.frame(c(pred_atuv[,i], pred_at_gj_un$ychains[,140+i], pred_at_mvco$ychains[,140+i]))
names(d_gg_bpso) <- c("pred")
#adding a column with the prediction types
d_gg_bpso$type <- c(rep("Univariate", length(pred_atuv[,i])),
                    rep("Unconditional Multivariate", length(pred_at_gj_un$ychains[,i])),
                    rep("Conditional Multivariate", length(pred_at_mvco$ychains[,i])))
d_gg_bpso$obs <- i
# Boxplot shows Boxes: 25% and 75 % Quantile; vertical line: median;
#lower whisker = smallest observation greater than or equal to lower hinge - 1.5 * IQR; 
#points: "remaining" outliers (die Erläuterungen komm in den Text unter der Abbildung)
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
#Unconditional and condtional predictions seem to behave similarly
