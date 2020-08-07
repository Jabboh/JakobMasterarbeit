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
d_res <- posterior_predict(fit_fin_cp, newdata = xdata, seed = 333)
#getting the predictions on the "probability scale, by taking the mean per 
#observation/column
pred <- colMeans(d_res)
#adding these predictions to data frame
ggd <- cbind(xdata, pred)
response_cp <- ggplot(data = ggd, aes(x = IA_500, y = pred, color = Mes)) +
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

#take the first 300 rows as predictions (deleting the last row)
pred_gj_un <- sim$sdList$yMu[1:300,1] 
dr_gj <- dr_gj[-nrow(dr_gj),]
#readding the Mes column
dr_gj$Mes <- ggd$Mes

ggd_gj_un <- cbind(dr_gj, pred_gj_un)
response_cp_gj <- ggplot(data = ggd_gj_un, aes(x = IA500, y = pred_gj_un, color = Mes)) +
  geom_point() + ylab("Predicted Probability") + xlab("IA_500 in Standard Units") +
  ggtitle ("Response Curve of Inundation Area for Culex Perexuus in GJAM") +
  guides(color=guide_legend(title="Month"))

#What to do now?:
#Überlegen wie man gjam mit rstanarm vergleichen kann?: (i) alles in einen Plot
# also Monate mit Farben und modeltype mit Zeichen, oder (ii) einen Monat
#aussuchen und dann die beiden Vergleichen
#2. Das Ganze noch mit conditional Predictions machen!
#3. Daraus eine Funktion bauen oder so, auf jeden Fall Code generischer machen!

