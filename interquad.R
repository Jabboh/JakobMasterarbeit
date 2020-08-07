###################WAIC Step Functions: 1. Dropping entire covariates and their
###################qassociated quadratic and interaction terms and 2. 
###################Dropping Interactions and quadratic terms

#Implementing Step 2 (Dropping interaction terms) and 3 (dropping quadratic
#terms).
interquad <- function(x, z){
  #define the base model to which the reduced models will be compared
  base <- x
  #defining the z vector without "mes", bc "mes" does not have interactions
  #nor quadratic terms
  z <- z[z != "Mes"]
  #creating a vector with the names of all the interactions
  int <- attr(base$terms, "term.labels")[grepl(":", attr(base$terms, "term.labels"))]
  #Create variable w_com that serves as the condition for the while-loop
  #After the first run, it is equal to the "best" model
  w_com <- matrix("start")
  rownames(w_com) <- "start"
  #initializing a while-loop that keeps dropping interactions as long as the 
  #less complex model is better. Second condition makes sure that there are 
  #actually interaction terms.
  while((rownames(w_com)[1] != "base") & (length(int) > 0)){
    #Create an empty list to store the model-objects of the different models
    st <- vector("list", length(int) + 1)
    #an empty list for all the WAICs
    waics <- st
    #add the baseline WAIC to which the more parsimonious models will 
    #be compared (It's always the "Most" complex model at the specific point of the 
    #model-selection procedure)
    waics[[length(waics)]] <- waic(base)
    #Calculate the different reduced models
    for (i in int){
      #index
      j <-match(i,int)
      #We run the model without the interaction terms of i and store the 
      #fitted model in st
      st[[j]] <- update(base, formula=drop.terms(base$terms, grep(i, attr(base$terms, "term.labels")),
                                                 keep.response=TRUE))
      print(attr(st[[j]]$terms, "term.labels"))
      #Calculating the WAIC of the model and storing it in waics
      waics[[j]] <- waic(st[[j]])
    }
    #a vector with all the model names (the ones without one covariate + the base)
    n <- c(int, "base")
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
    #deleting the dropped variable in the names vector (Es sind ja immer zwei,
    #weil ja zwei Variablen mit einem Interaction term assoziert sind)
    int <- int[int != rownames(w_com)[1]] 
  }
  #Now doing more or less the same with the quadratic terms
  #resetting w_com variable
  w_com <- matrix("start")
  rownames(w_com) <- "start"
  #initializing a while-loop that keeps dropping interactions as long as the 
  #one less complex model is better
  while(rownames(w_com)[1] != "base"){
    #Create an empty list to store the model-objects of the different models
    st <- vector("list", length(z) + 1)
    #an empty list for all the WAICs
    waics <- st
    #add the baseline WAIC to which the more parsimonious models will 
    #be compared (It's always the "Most" complex model at the specific point of the 
    #model-selection procedure)
    waics[[length(waics)]] <- waic(base)
    #Calculate the different reduced models
    for (i in z){
      #index
      j <-match(i,z)
      #We run the model without the terms associated with z (linear, quadratic
      #and interaction terms and store the fitted model in st)
      st[[j]] <- update(base, formula=drop.terms(base$terms, grep(paste0("I(", i), attr(base$terms, "term.labels"), fixed = T),
                                                 keep.response=TRUE))
      print(attr(st[[j]]$terms, "term.labels"))
      #Calculating the WAIC of the model and storing it in waics
      waics[[j]] <- waic(st[[j]])
    }
    #a vector with all the model names (the ones without one covariate + the base)
    n <- c(z, "base")
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
  base#for some reason this does not work
}