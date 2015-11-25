### Functions solving calibration with minimum bounds (for bounded distances).
### Are all private, used by the main "calibration" function

minBoundsCalib <- function(Xs, d, total, q=rep(1,length(d)),
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
