### Functions solving calibration with minimum bounds (for bounded distances).
### Are all private, used by the main "calibration" function

solveMinBoundsCalib <- function(Xs, d, total, q=NULL,
                           maxIter=500, calibTolerance=1e-06, description=TRUE) {

  if (!requireNamespace("lpSolve", quietly = TRUE)) {
      stop("Package lpSolve needed for this function to work. Please install it.",
            call. = FALSE)
  }

  n <- length(d)
  oneN <- rep(1,n)
  zeroN <- rep(0,n)
  IdentityN <- diag(1,n)
  B <- t(diag(d) %*% Xs)

  a <- c(1,rep(0,n))

  A1 <- rbind( cbind(-oneN,IdentityN) , cbind(oneN,IdentityN) )
  b1 <- c(oneN,oneN)

  A3 <- matrix(cbind(rep(0,nrow(B)+1),rbind(B,zeroN)), ncol= n+1, nrow= nrow(B)+1)
  b3 <- c(total,0)

  Amat <- rbind(A1,A3)
  bvec <- c(b1,b3)
  const <- c(rep("<=",n), rep(">=",n), rep("=", length(b3)))

#   simplexSolution <- simplex(a, A1=A1, b1=b1, A3=A3, b3=b3)
#   simplexSolution <- solveLP(a, bvec=bvec, Amat=Amat, const.dir=const, lpSolve=T,
#                               maxiter=10000, tol=1e-2, verbose=1)
  simplexSolution <- lpSolve::lp(direction="min", objective.in = a, const.mat = Amat,
                          const.dir = const, const.rhs = bvec)

  xSol <- simplexSolution$solution
  minBounds <- xSol[1]
  gSol <- xSol[2:(n+1)]
  wSol <- gSol * d

  if(description) {
    writeLines("Solution found for calibration on minimal bounds:")
    writeLines(paste("L =",min(xSol)))
    writeLines(paste("U =",max(xSol)))
  }

  ## TODO : what if there is no solution ?
  ## Test avec données ex2:
  # 1.663102
  # 0.3368984

  return(gSol)

}

minBoundsCalib <- function(Xs, d, total, q=NULL,
                           maxIter=500, calibTolerance=1e-06, description=TRUE,
                           precisionBounds=1e-4, forceSimplex=FALSE) {


  if(forceSimplex || nrow(Xs) <= 1000) {
    
    gSol <- solveMinBoundsCalib(Xs, d, total, q,
                                maxIter, calibTolerance, description)
    
    Lmax <- min(gSol)
    Umin <- max(gSol)
    
  } else {
    
    Lmax <- 1.0
    Umin <- 1.0
    
  }


  digitsPrec <- abs(log(precisionBounds,10))

#   print("Bornes test :")
  Ltest1 <- round(Lmax - 5*10**(-digitsPrec-1),digitsPrec)
  Utest1 <- round(Umin + 5*10**(-digitsPrec-1),digitsPrec)

  Ltest <- Ltest1
  Utest <- Utest1

#   print(Ltest)
#   print(Utest)

  gFinal <- calib(Xs=Xs, d=d, total=total, method="logit", bounds = c(Ltest,Utest),
                  maxIter=maxIter, calibTolerance=calibTolerance)

  ## If no convergence, bisection to find the true min Bounds
  if(is.null(gFinal)) {

    gTestMin <- calib(Xs=Xs, d=d, total=total, method="raking", maxIter=maxIter, calibTolerance=calibTolerance)
    LtestMin <- min(gTestMin)
    LtestMax <- max(gTestMin)

    return(bisectMinBounds(c(LtestMin,LtestMax),c(Ltest1,Utest1),gTestMin,
                                       Xs,d,total, method,maxIter, calibTolerance, precisionBounds, description))

    Ltest <- Ltest1
    Utest <- Utest1

  }

  return(gFinal)

}

bisectMinBounds <- function(convergentBounds,minBounds,gFinalSauv,
                            Xs,d,total, method,maxIter, calibTolerance, precisionBounds,
                            description=TRUE) {

  newBounds <- (minBounds + convergentBounds) / 2
  Ltest <- newBounds[1]
  Utest <- newBounds[2]

  if(description) {
    writeLines(paste("------ Bisection search : ",newBounds[1],";",newBounds[2]))
  }

  gFinal <- calib(Xs=Xs, d=d, total=total, method="logit", bounds = c(Ltest,Utest),
                   maxIter=maxIter, calibTolerance=calibTolerance)

  if(is.null(gFinal)) {
    return( bisectMinBounds(convergentBounds,c(Ltest,Utest),gFinalSauv,
                            Xs,d,total, method,maxIter, calibTolerance, precisionBounds, description) )
  } else {

    if( all(abs(gFinal - gFinalSauv) <= rep(precisionBounds, length(gFinal))) ) {
      return(gFinal)
    } else {
      return( bisectMinBounds(c(Ltest,Utest),minBounds,gFinal,
                              Xs,d,total, method,maxIter, calibTolerance, precisionBounds, description) )
    }

  }

}
