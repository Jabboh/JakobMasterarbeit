interactions <- function(x, y){
  #define the base model to which the reduced models will be compared
  base <- x
  #defining the y vector with "mes", bc "mes" does not have interactions
  y <- y[y != "Mes"]
  #Create variable w_com that serves as the condition for the while-loop
  #After the first run, it is equal to the "best" model
  
  #resetting w_com variable
  w_com <- matrix("start")
  rownames(w_com) <- "start"
  #initializing a while-loop that keeps dropping interactions as long as the 
  #one less complex model is better
  while((rownames(w_com)[1] != "base") & (grepl(":", attr(base$terms, "term.labels")) %>% any)){
    #Create an empty list to store the model-objects of the different models
    st <- vector("list", length(y) + 1)
    #an empty list for all the WAICs
    waics <- st
    #add the baseline WAIC to which the more parsimonious models will 
    #be compared (It's always the "Most" complex model at the specific point of the 
    #model-selection procedure)
    waics[[length(waics)]] <- waic(base)
    #Calculate the different reduced models
    for (i in y){
      #index
      j <-match(i,y)
      #We run the model without the terms associated with y (linear, quadratic
      #and interaction terms and store the fitted model in st)
      st[[j]] <- update(base, formula=drop.terms(base$terms, grep(paste0(i, ".*:|:.*", i), attr(base$terms, "term.labels")),
                                                 keep.response=TRUE))
      print(attr(st[[j]]$terms, "term.labels"))
      #Calculating the WAIC of the model and storing it in waics
      waics[[j]] <- waic(st[[j]])
    }
    #a vector with all the model names (the ones without one covariate + the base)
    n <- c(y, "base")
    #naming the two storage lists accordingly
    names(st) <- n
    names(waics) <- n
    #Compare the different WAICs, The first row contains the preferred model
    #(the one with the lowest WAIC)
    w_com <- loo_compare(waics)
    print(w_com)
    #Printing a message to the screen indicating that variable xxx will be thrown out
    print(paste0 ("The interaction terms of variable ", rownames(w_com)[1], " will be dropped"))
    
    #storing the model without the removed variable as the new base model
    base <- st[[which(names(st) == rownames(w_com)[1])]]
  }
  #Return the final model
  return(base)#for some reason this does not work
}

arg <- interactions(x,y)

#Now, I check whether the quadratic terms, "improve the model"
quadratic <- function(x, y){
  #define the base model to which the reduced models will be compared
  base <- x
  #resetting w_com variable
  w_com <- matrix("start")
  rownames(w_com) <- "start"
  #initializing a while-loop that keeps dropping interactions as long as the 
  #one less complex model is better
  while(rownames(w_com)[1] != "base"){
    #Create an empty list to store the model-objects of the different models
    st <- vector("list", length(y) + 1)
    #an empty list for all the WAICs
    waics <- st
    #add the baseline WAIC to which the more parsimonious models will 
    #be compared (It's always the "Most" complex model at the specific point of the 
    #model-selection procedure)
    waics[[length(waics)]] <- waic(base)
    #Calculate the different reduced models
    for (i in y){
      #index
      j <-match(i,y)
      #We run the model without the terms associated with y (linear, quadratic
      #and interaction terms and store the fitted model in st)
      st[[j]] <- update(base, formula=drop.terms(base$terms, grep(paste0("I(", i), attr(base$terms, "term.labels"), fixed = T),
                                                 keep.response=TRUE))
      print(attr(st[[j]]$terms, "term.labels"))
      #Calculating the WAIC of the model and storing it in waics
      waics[[j]] <- waic(st[[j]])
    }
    #a vector with all the model names (the ones without one covariate + the base)
    n <- c(y, "base")
    #naming the two storage lists accordingly
    names(st) <- n
    names(waics) <- n
    #Compare the different WAICs, The first row contains the preferred model
    #(the one with the lowest WAIC)
    w_com <- loo_compare(waics)
    print(w_com)
    #Printing a message to the screen indicating that variable xxx will be thrown out
    print(paste0 ("The quadratic term of variable ", rownames(w_com)[1], " will be dropped"))
    
    #storing the model without the removed variable as the new base model
    base <- st[[which(names(st) == rownames(w_com)[1])]]
  }
  #Return the final model
  return(base)#for some reason this does not work
}
