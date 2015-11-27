#### All functions in this file are private methods used by the
#### "calibration" function

# TODO : each inverseDistance function comes with an "updateParameters" function

# Fonction Wrapper
calib <- function(Xs, d, total, q=NULL, method=NULL, bounds = NULL,
                  alpha = NULL,
                  maxIter=500, calibTolerance=1e-06) {

  if(!is.null(method)) {
  switch(method,
          linear={
            inverseDistance <- inverseDistanceLinear
            params <- NULL
            updateParameters <- identity
          },
          raking={
            inverseDistance <- inverseDistanceRaking
            params <- NULL
            updateParameters <- identity
          },
          logit={
            inverseDistance <- inverseDistanceLogit
            params <- bounds
            updateParameters <- identity
          },
           truncated={
             inverseDistance <- inverseDistanceTruncated
             params <- bounds
             updateParameters <- identity
           },
          curlingHat={
            inverseDistance <- distanceCurlingHat
            # TODO : check params in list are correctly entered
            params <- c(0.5,1.5) # For tests only
            updateParameters <- updateParametersCurlingHat
          },
          {
            print('By default, raking method selected')
            params <- NULL
            inverseDistance <- inverseDistanceRaking
            updateParameters <- identity
          }
  )
  } else {
    print('By default, raking method selected')
    params <- NULL
    inverseDistance <- inverseDistanceRaking
    updateParameters <- identity
  }
  # TODO : additional checks ??

  return(calibAlgorithm(Xs, d, total, q, inverseDistance,
                        updateParameters, params, maxIter, calibTolerance))

}

calibAlgorithm <- function(Xs, d, total, q=NULL,
                            inverseDistance, updateParameters, params,
                            maxIter=500, calibTolerance=1e-06) {

  if(is.null(q)) {
    q <- rep(1,length(d))
  }

  g <- NULL
  toleranceGInv = .Machine$double.eps # Tolerance when we compute ginv

  ## Linear optimization algorithm
  lambda = as.matrix(rep(0, ncol(Xs)))
  wTemp = as.vector(d * inverseDistance(Xs %*% lambda * q, params))

  cont <- TRUE
  l <- 1

  while (cont) {

    phi = t(Xs) %*% wTemp - total
    T1 = t(Xs * wTemp)
    phiprim = T1 %*% Xs
    
    wTemp <- NA
    
    try({
      lambda = lambda - ginv(phiprim, tol = toleranceGInv) %*% phi
      wTemp = as.vector(d * inverseDistance(Xs %*% lambda * q, params))
      }
    )

    if (any(is.na(wTemp)) | any(is.infinite(wTemp))) {
      warning("No convergence")
      return(NULL)
    }

    tHat = t(Xs) %*% wTemp
    if (max(abs(tHat - total)/total) < calibTolerance) {
      cont <- FALSE
    }

    if(l >= maxIter) {
      cont <- FALSE
      warning(paste("No convergence in ", maxIter, " iterations."))
      return(NULL)
    }

    l <- l+1
    # update Parameters before going back into loop
    updateParameters(params)

  }

  ## Return solution if found
  if(l <= maxIter) {
    g = wTemp/d
    return(g)
  } else {
    warning(paste("No convergence in ", maxIter, " iterations."))
    return(NULL)
  }

  return(g)
}


# TODO : add qk vectors

# Explain why params (only use S3 and not S4 methods)
inverseDistanceLinear <- function(x, params=NULL) {
  return(1+x)
}

inverseDistanceRaking <- function(x, params=NULL) {
  return(exp(x))
}

# Params are bounds
inverseDistanceLogit <- function(x, bounds) {

  if(length(bounds) != 2) {
    stop("Must enter LO and UP bounds in a vector.
          Example : bounds=c(0.5,1.5) for LO = 0.5 and UP=1.5")
  }

  L = bounds[1]
  U = bounds[2]
  A = (U - L) / ((1-L)*(U-1))
  distance = ( L*(U-1) + U*(1-L)*exp(A*x) ) / ( (U-1) + (1-L)*exp(A*x))

  return(distance)
}

# TODO : truncated method ?
# TODO : hyperbolic sine

# TODO : shaped -> "param" function which updates parameters is
# contained within distance function
