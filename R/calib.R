
# TODO : Replace bounds by more generic "distanceParams"
# TODO : remove description=FALSE everywhere
# Params should be a list containing pointers to useful parameters

# TODO : each distance function comes with an "updateParameters" function


# TODO : fonction pour faire le calage avec des paramÃ¨tres et une fonction de distance,
# et fonction "wrapper", qui sera celle appelÃ©e par l'utilisateur.

# Fonction Wrapper
calib <- function(Xs, d, total, q=rep(1,length(d)), method=NULL, bounds = NULL,
                  alpha = NULL,
                  maxIter=500, calibTolerance=1e-06) {

  if(!is.null(method)) {
  switch(method,
          linear={
            distance <- pseudoDistanceLinear
            params <- NULL
            updateParameters <- identity
          },
          raking={
            distance <- pseudoDistanceRaking
            params <- NULL
            updateParameters <- identity
          },
          logit={
            distance <- pseudoDistanceLogit
            # TODO : check params in list are correctly entered
            params <- bounds
            updateParameters <- identity
          },
          curlingHat={
            distance <- distanceCurlingHat
            # TODO : check params in list are correctly entered
            params <- c(0.5,1.5) # For tests only
            updateParameters <- updateParametersCurlingHat
          },
          {
            print('By default, raking method selected')
            params <- NULL
            distance <- pseudoDistanceRaking
            updateParameters <- identity
          }
  )
  } else {
    print('By default, raking method selected')
    params <- NULL
    distance <- pseudoDistanceRaking
    updateParameters <- identity
  }
  # TODO : additional checks ??

  return(calibAlgorithm(Xs, d, total, q, distance,
                        updateParameters, params, maxIter, calibTolerance))

}

calibAlgorithm <- function (Xs, d, total, q=rep(1,length(d)),
                            distance, updateParameters, params,
                            maxIter=500, calibTolerance=1e-06) {

  g <- NULL
  toleranceGInv = .Machine$double.eps # Tolerance when we compute ginv

  ## Linear optimization algorithm
  lambda = as.matrix(rep(0, ncol(Xs)))
  wTemp = as.vector(d * distance(Xs %*% lambda * q, params))

  cont <- TRUE
  l <- 1

  while (cont) {

    phi = t(Xs) %*% wTemp - total
    T1 = t(Xs * wTemp)
    phiprim = T1 %*% Xs
    lambda = lambda - ginv(phiprim, tol = toleranceGInv) %*% phi
    wTemp = as.vector(d * distance(Xs %*% lambda * q, params))


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


# TODO : explain clearly that these are 
# not real distances, but inverse of distances
# TODO : add qk vectors
# TODO : separate file for distances ?

# Explain why params (only use S3 and not S4 methods)
pseudoDistanceLinear <- function(x, params=NULL) {
  return(1+x)
}

pseudoDistanceRaking <- function(x, params=NULL) {
  return(exp(x))
}

# Params are bounds
pseudoDistanceLogit <- function(x, bounds) {

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

# TODO : truncated linear
# TODO : hyperbolic sine

# TODO : shaped -> "param" function which updates parameters is
# contained within distance function
