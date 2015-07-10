
# Computes the Horvitz-Thompson estimator for
# the total of @param var column (given @param weights)
HTtotal <- function(var, weights) {
  
  if(!is.numeric(var)) {
    var <- as.numeric(var)
  }
  
  if(!is.numeric(weights)) {
    weights <- as.numeric(weights)
  }
  
  
  return(var %*% weights)
}

HTmean <- function(var, weights, popTot=NULL) {
  
  # Population total defaults to sum of weights
  if(is.null(popTot)) {
    
    if(is.numeric(weights)) {
      popTot <- sum(weights)
    } else {
      popTot <- sum(as.numeric(weights))
      warning("weights column is not numeric. Automatic conversion done ; 
              could lead to bugs")
    }
    
  }
  
  return( HTtotal(var, weights) / popTot )
}