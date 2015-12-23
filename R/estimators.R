# copyright (C) 2015 A.Rebecq

#' Computes the Horvitz-Thompson estimator for the total of a column
#' @param var column of variable of interest
#' @param weights column of weights matching the variable of interest
#' @return Estimated total
#' @export
HTtotal <- function(var, weights) {

  if(!is.numeric(var)) {
    var <- as.numeric(var)
  }

  if(!is.numeric(weights)) {
    weights <- as.numeric(weights)
  }


  return(var %*% weights)
}

#' Computes the Horvitz-Thompson estimator for the mean of a column
#' @param var column of variable of interest
#' @param weights column of weights matching the variable of interest
#' @return Estimated mean
#' @export
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
