###Response Curve für IA_500
#a sequence of IA_500-values from min to max with 1000 steps
ia <- seq(min(df$IA_500), max(df$IA_500), length.out = 50)
#names of covariates
nc <- c("IA_500", "NDVIBEF_2000", "Mes")
#the mean of the covariates is 0, bc I normalized them
#predict the probability of presence for the different ia-values, holding the
#other covariates at their mean (0) and fixing mes = abril
#creating  x-data dataframe with all the avareges (which are zero)
data <- as.data.frame(matrix(0, ncol = length(nc), nrow = length(ia)))
names(data) <- nc  
#replace IA_500 0s with the sequence from above
data$IA_500 <-ia
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
d_res <- posterior_linpred(fit_fin_cp, transform = T, newdata = xdata, seed = 333)
#getting the predictions on the "probability scale, by taking the mean per 
#observation/column
univariate <- colMeans(d_res)
#adding these predictions to data frame
ggd <- cbind(xdata, univariate)
response_cp <- ggplot(data = ggd, aes(x = IA_500, y = univariate, color = Mes)) +
  geom_point() + ylab("Predicted Probability") + xlab("IA_500 in Standard Units") +
  ggtitle ("Response Curve of Inundation Area for Culex Perexuus in Univariate Model") +
  guides(color=guide_legend(title="Month"))

##Doing the same thing for gjam

#for unconditional in-sample predictions
#Converting the Mes-Variable to a bunch of dummies

#create a storage matrix
st <- data.frame(matrix(, nrow=nrow(xdata), ncol = 6))
#assign the factor names to the column names
names(st) <- levels(df$Mes)
#loop over all the levels
for (i in levels(xdata$Mes)){
  #loop over all the observations
  for (j in 1:nrow(xdata)){
    st[j, i] <- ifelse(xdata[j, "Mes"] == i, 1, 0) #fill the dummy data frame
  }
}

#combinig our original data with dummy data frame
dr_gj<- cbind(xdata, st)
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

ggd_gj_un <- cbind(dr_gj, pred_gj_un)
response_cp_gj <- ggplot(data = ggd_gj_un, aes(x = IA500, y = pred_gj_un, color = Mes)) +
  geom_point() + ylab("Predicted Probability") + xlab("IA_500 in Standard Units") +
  ggtitle ("Response Curve of Inundation Area for Culex Perexuus in GJAM") +
  guides(color=guide_legend(title="Month"))

#Plot gjam and rstanarm response curves in the same plot

#prepare the dataframe to feed into ggplot (just add predictions of gjam to dataframe of rstanarm)
gg_cp_uv_un <- ggd
gg_cp_uv_un$multivariate <- ggd_gj_un$pred_gj_un

#Do the ggplot 

response_cp_uv_un <- ggplot(data = gg_cp_uv_un, aes(x = IA_500, color = Mes)) +
  geom_point(aes(y = univariate, shape = "univariate")) + 
  geom_point(aes(y = multivariate, shape = "multivariate")) + ylab("Predicted Probability") + 
  xlab("IA_500 in Standard Units") +
  ggtitle ("Univariate vs. Multivariate Response Curves of Inundation Area for Culex Perexuus") +
  guides(color=guide_legend(title="Month"), shape = guide_legend(title="Model"))

#Why are the response curves different? Does that contradict our hypothesis that unconditional
#predictions dont differ from univariate predictions?
#What to do now?:
#Überlegen wie man gjam mit rstanarm vergleichen kann?: (i) alles in einen Plot
# also Monate mit Farben und modeltype mit Zeichen, oder (ii) einen Monat
#aussuchen und dann die beiden Vergleichen
#2. Das Ganze noch mit conditional Predictions machen!
#3. Daraus eine Funktion bauen oder so, auf jeden Fall Code generischer machen!

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
ydata <- tibble(c(rep(0, 300), rep(1, 300)))
names(ydata) <- "Anatropre"
ydata <- rbind(ydata,1)

#define the modelling settings for posterior simulations with the conditional species as present
newdata <- list(xdata = dr_gj_con, ydataCond = ydata, nsim = 4000)

#calculating the in-sample predictions (simulations)
sim <- gjamPredict(output = joint_fin, newdata = newdata)

#take the first 300 rows as predictions (deleting the last row); Predictions are the means of the
#ychains
pred_gj_con <- sim$sdList$yMu[1:600,1] 

#adding the conditional presence of Anotrophes
dr_gj_con$cond <- factor(ydata$Anatropre, labels = c("absent", "present"))
#deleting the last row of our data (remember: we added only bc this way the predict function runs)
dr_gj_con <- dr_gj_con[-601,]
#readding the Mes column
dr_gj_con$Mes <- ggd$Mes

ggd_gj_con <- cbind(dr_gj_con, pred_gj_con)
response_cp_gj_con <- ggplot(data = ggd_gj_con, aes(x = IA500, y = pred_gj_con, color = Mes, shape = cond)) +
  geom_point() + ylab("Predicted Probability") + xlab("IA_500 in Standard Units") +
  ggtitle ("Response Curve of Inundation Area for Culex Perexuus in GJAM \n Conditional on Presence of Anophles Atroparvus") +
  guides(color=guide_legend(title="Month"), shape = guide_legend(title = "Anopheles atroparvus"))

#strange results, ordering of the effect of month on pred, changes depending on whether AT is present
#or not. Does that make sense? Is that even possible? >> I actually do not think so (BUt only september
#is different!)

###Doing a plot of gjam conditional predictions and univariate prediction (We just take any month)

#taking June


#Preparing the input data
#number of "Observations
nu <- nrow(xdata[xdata["Mes"] == "Junio",])
d <- rbind(xdata[xdata["Mes"] == "Junio",], xdata[xdata["Mes"] == "Junio",], xdata[xdata["Mes"] == "Junio",])
#creating one column that has the presence-absence of at for the gjam models and for the univariate
#model the variable takes on another factor level (2)
d$mode <- factor(c(rep(0, nu), rep(1, nu), rep(2, nu)), labels = c("GJAM with absent Anopheles", 
                                                                   "GJAM with present Anopheles",
                                                                   "Univariate Model"))
#adding the predictions
d$pred <- c(ggd_gj_con$pred_gj_con[ggd_gj_con$Mes == "Junio"], ggd$univariate[ggd$Mes == "Junio"])

#Doing the ggplot
response_cp_gjuni_con <- ggplot(data = d, aes(x = IA_500, y = pred, shape = mode)) +
  geom_point() + ylab("Predicted Probability") + xlab("IA_500 in Standard Units") +
  ggtitle ("Response Curve of Inundation Area for Culex Perexuus in GJAM \n Conditional on Presence of Anophles Atroparvus vs. Univariate Model for Month June") +
  guides(shape = guide_legend(title = "Anopheles atroparvus"))

#man könnte eine Reihe solcher Plots für jeden Monat machen

#Plots, die ich haben möchte: alle von oben (beim letzten eine Reihe über alle Monate) + nochmal
#einen Plot die unconditional gjam response curves on their own

#Try to make a function out of it (But, if it takes too long, then do them by hand separately)
#>> I think that's too cumbersome, rather make your code easier and more tracable

