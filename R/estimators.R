# copyright (C) 2014-2016 A.Rebecq

#' Weighted estimator for total
#' @description 
#' Computes the weighted estimator for the total of a column. Alias for
#' \code{\link{weightedTotal}}
#' @param var column of variable of interest
#' @param weights column of weights matching the variable of interest
#' @return Estimated total
#' @seealso \code{\link{weightedTotal}}
#' @export
HTtotal <- function(var, weights) {
  return(weightedTotal(var, weights))
}

#' Weighted estimator for the mean
#' @description 
#' Computes the weighted estimator for the mean of a column. Alias for
#' \code{\link{weightedMean}}
#' @param var column of variable of interest
#' @param weights column of weights matching the variable of interest
#' @return Estimated mean
#' @seealso \code{\link{weightedMean}}
#' @export
HTmean <- function(var, weights, popTot=NULL) {
 return(weightedMean(var, weights, popTot)) 
}

#' Weighted estimator for total
#' @description 
#' Computes the weighted estimator for the total of a column
#' @param var column of variable of interest
#' @param weights column of weights matching the variable of interest
#' @return Estimated total
#' @seealso \code{\link{HTtotal}}
#' @export
weightedTotal <- function(var, weights) {

  if(!is.numeric(var)) {
    var <- as.numeric(var)
  }

  if(!is.numeric(weights)) {
    weights <- as.numeric(weights)
  }


  return(var %*% weights)
}

#' Weighted estimator for the mean
#' @description 
#' Computes the weighted estimator for the mean of a column
#' @param var column of variable of interest
#' @param weights column of weights matching the variable of interest
#' @return Estimated mean
#' @seealso \code{\link{HTmean}}
#' @export
weightedMean <- function(var, weights, popTot=NULL) {

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
