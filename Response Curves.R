##For Culex Perexuus

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
  ggtitle ("Response Curve of Inundation Area for Culex Perexuus in Univariate Model") +
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
  ggtitle ("Response Curve of Inundation Area for Culex Perexuus in GJAM") +
  guides(color=guide_legend(title="Month"))

##Plot gjam and rstanarm response curves in the same plot

#prepare the dataframe to feed into ggplot (just add predictions of gjam to dataframe of rstanarm)

ggd$multivariate <- ggd_gj_un$pred_gj_un

#Do the ggplot 

response_cp_uv_un <- ggplot(data = ggd, aes(x = IA_500, color = Mes)) +
  geom_point(aes(y = univariate, shape = "univariate")) + 
  geom_point(aes(y = multivariate, shape = "multivariate")) + ylab("Predicted Probability") + 
  xlab("IA_500 in Standard Units") +
  ggtitle ("Univariate vs. Multivariate Response Curves of Inundation Area for Culex Perexuus") +
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
  ggtitle ("Response Curve of Inundation Area for Culex Perexuus in GJAM \n Conditional on Presence of Anophles Atroparvus") +
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
    ggtitle (paste0("Response Curve of Inundation Area for Culex Perexuus in GJAM \n Conditional on Presence of Anophles Atroparvus vs. Univariate Model for Month ", i)) +
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
  ggtitle ("Response Curve of NDVI Before for Culex Perexuus in Univariate Model") +
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
  ggtitle ("Response Curve of NDVI Before for Culex Perexuus in GJAM") +
  guides(color=guide_legend(title="Month"))

##Plot gjam and rstanarm response curves in the same plot

#prepare the dataframe to feed into ggplot (just add predictions of gjam to dataframe of rstanarm)

ggd$multivariate <- ggd_gj_un$pred_gj_un

#Do the ggplot 

response_cp_uv_un <- ggplot(data = ggd, aes(x = NDVIBEF_2000, color = Mes)) +
  geom_point(aes(y = univariate, shape = "univariate")) + 
  geom_point(aes(y = multivariate, shape = "multivariate")) + ylab("Predicted Probability") + 
  xlab("NDVIBEF_2000 in Standard Units") +
  ggtitle ("Univariate vs. Multivariate Response Curves of NDVI Before for Culex Perexuus") +
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
  ggtitle ("Response Curve of NDVI Before for Culex Perexuus in GJAM \n Conditional on Presence of Anophles Atroparvus") +
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
    ggtitle (paste0("Response Curve of NDVI Before for Culex Perexuus in GJAM \n Conditional on Presence of Anophles Atroparvus vs. Univariate Model for Month ", i)) +
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
  ggtitle ("Univariate vs. Multivariate Response Curves of NDVI Before for Anopheles Troparvus") +
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