#make a function that converts a factor variable "x" into levels-1-dummy variables that gjam can 
#handle in its modell fitting and predictions. y is the dataframe of data to which the dummies should
#added

gjamDummy <- function(x, y){
  #create a storage matrix
  st <- data.frame(matrix(, nrow=length(x), ncol = x %>% levels %>% length))
  names(st) <- x %>% levels
  #loop over all the factor levels and store them as dummies in separate columns of st
  for (i in levels(x)){
    #loop over all the observations
    for (j in 1:length(x)){
      st[j, i] <- ifelse(x[j] == i, 1, 0) #fill the dummy data frame
    }
  }
  #combinig our original data with dummy data frame
  y<- cbind(y, st)
  return (y)
}
