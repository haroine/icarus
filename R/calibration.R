
# This function executes easy calibration with just data and matrix of margins
# TODO : add calib options
# - scale : if TRUE, weights ratio are centered on 1 (just like "ECHELLE=0" in Calmar2)
# - description : print summary of before / after weight ratios and comparison
# of estimators before / after calibration (TRUE by default)
# To only calibrate on popTotal, set marginMatrix=NULL.
# TODO : remove useless parameter "description"
# Careful, contrary to old rule, returns calibrated weights and not
# ratio between initial and calibrated weights
# Careful with lambda, triggers restrictedSearch, meaning searchLambda could never converge
#########
#' Calibration on margins
#' @description
#' Performs calibration on margins with several methods and customizable parameters
#' @param data The dataframe containing the survey data
#' @param marginMatrix The matrix giving the margins for each column variable included
#' in the calibration problem
#' @param colWeights The name of the column containing the initial weights in the survey
#' dataframe
#' @param colCalibratedWeights The name of the column of final calibrated weights
#' @param method The method used to calibrate. Can be "linear", "raking", "logit", "truncated"
#' @param maxIter The maximum number of iterations before stopping
#' @param description If TRUE, output stats about the calibration process as well as the
#' graph of the density of the ratio calibrated weights / initial weights
#' @param bounds Two-element vector containing the lower and upper bounds for bounded methods
#' ("truncated" and "logit")
#' @param costs The penalized calibration method will be used, using costs defined by this
#' vector. Must match the number of rows of marginMatrix. Negative of non-finite costs are given
#' a default that can be adjusted through parameter "infinity"
#' @param popTotal Precise the total population if margins are defined by relative value in
#' marginMatrix (percentages)
#' @param scale If TRUE, stats (including bounds) on ratio calibrated weights / initial weights are
#' done on a vector multiplied by the weighted non-response ratio (ratio population total /
#' total of initial weights). Has same behavior as "ECHELLE=0" in Calmar.
#' @param check performs a few check about the dataframe. TRUE by default
#' @param infinity Only used in the penalized calibration. Use this to tweak the numeric value
#' of an infinite cost.
#' @param uCostPenalized Unary cost by which every cost is "costs" column is multiplied
#' @param lambda The initial ridge lambda used in penalized calibration. By default, the initial
#' lambda is automatically chosen by the algorithm, but you can speed up the search for the optimum
#' if you already know a lambda close to the lambda_opt corresponding to the gap you set. Be careful,
#' the search zone is reduced when a lambda is set by the user, so the program may not converge
#' if the lambda set is too far from the lambda_opt.
#' @param gap Only useful for penalized calibration. Sets the maximum gap between max and min
#' calibrated weights / initial weights ratio (and thus is similar to the "bounds"
#' parameter used in regular calibration)
#' @param precisionBounds Only used for calibration on minimum bounds. Desired precision
#' for lower and upper reweighting factor, both bounds being as close to 1 as possible
#' @param forceSimplex Only used for calibration on tight bounds.Bisection algorithm is used
#' for matrices whose size exceed 1e8. forceSimplex = TRUE forces the use of the simplex algorithm
#' whatever the size of the problem (you might want to set this parameter to TRUE if you
#' have a large memory size)
#' @param exportDistributionImage File name to which the density plot shown when
#' description is TRUE is exported. Requires package "ggplot2"
#' @param exportDistributionTable File name to which the distribution table of before/after
#' weights shown when description is TRUE is exported.
#' Requires package "ggplot2". Requires package "xtable"
#'
#' @return column containing the final calibrated weights
#'
#' @export
calibration = function(data, marginMatrix, colWeights = "POIDS", colCalibratedWeights="POIDS_CALES", method="linear",
                       maxIter=2500, description=TRUE, bounds=NULL, costs=NULL, popTotal=NULL, scale=NULL, check=TRUE
                       , infinity=1e7, uCostPenalized=1e2, lambda=NULL, gap=NULL, precisionBounds=1e-4, forceSimplex=FALSE
                       , exportDistributionImage=NULL, exportDistributionTable=NULL) {
  
  # By default, scale is TRUE when popTotal is not NULL, false otherwise
  if(is.null(popTotal)) {
    scale <- FALSE
  } else {
    scale <- TRUE
  }
  
  ## Clean costs for penalized calibration
  #   if(!is.null(costs)) {
  #     costs <- cleanCosts(costs, infinity, uCostPenalized)
  #   }
  
  if(check) {
    
    # Check NAs on calibration variables
    matrixTestNA = missingValuesMargins(data, marginMatrix)
    testNA = as.numeric(matrixTestNA[,2])
    if(sum(testNA) > 0) {
      print(matrixTestNA)
      stop("NAs found in calibration variables")
    }
    
    # check if number of modalities in calibration variables matches marginMatrix
    if(!checkNumberMargins(data, marginMatrix)) stop("Error in number of modalities.")
    
  }
  
  ## TODO : add parameter pct
  marginCreation <- createFormattedMargins(data, marginMatrix, popTotal)
  matrixCal = marginCreation[[2]]
  formattedMargins = marginCreation[[1]]
  
  # Same rule as in "Calmar" for SAS : if scale is TRUE,
  # calibration is done on weights adjusted for nonresponse
  # (uniform adjustment)
  weights <- as.numeric(data.matrix(data[colWeights]))
  
  ## For tests
  #   Xs_glob <<- matrixCal
  #   total_glob <<- formattedMargins
  #   d_glob <<- weights
  
  if(scale) {
    if(is.null(popTotal)) {
      stop("When scale is TRUE, popTotal cannot be NULL")
    }
    weights <- weights*(popTotal / sum(data.matrix(data[colWeights])) )
    
  }
  
  # TODO : add catching of "no convergence" exception
  if(is.null(costs)) {
    
    g <- NULL
    
    if( (is.numeric(bounds)) || (method != "min") ) {
      g <- calib(Xs=matrixCal, d=weights, total=formattedMargins, method=method, bounds=bounds, maxIter=maxIter)
    } else {
      if( (bounds == "min") || (method == "min")) {
        g <- minBoundsCalib(Xs=matrixCal, d=weights, total=formattedMargins
                            , q=rep(1,length(d)), maxIter=maxIter, description=description, precisionBounds=precisionBounds, forceSimplex=forceSimplex)
      }
    }
    
    data[colCalibratedWeights] = g*weights
    
  } else {
    
    # Forbid popTotal null when gap is selected
    if(!is.null(gap) && is.null(popTotal)) {
      warning("popTotal NULL when gap is selected is a risky setting !")
    }
    
    ## Format costs
    costsFormatted <- formatCosts(costs, marginMatrix, popTotal)
    
    wCal = penalizedCalib(Xs=matrixCal, d=weights, total=formattedMargins, method=method
                          , bounds=bounds, costs=costsFormatted, infinity=infinity, uCostPenalized=uCostPenalized
                          , maxIter=maxIter, lambda=lambda, gap=gap)
    data[colCalibratedWeights] = data.matrix(wCal)
    g = wCal / weights
  }
  
  if(description) {
    writeLines("")
    writeLines("################### Summary of before/after weight ratios ###################")
  }
  # popTotalComp is popTotal computed from sum of calibrated weights
  popTotalComp <- sum(data[colCalibratedWeights])
  
  weightsRatio = g
  
  if(description) {
    
    writeLines(paste("Calibration method : ",method, sep=""))
    
    if(! (method %in% c("linear","raking")) && ! is.null(bounds) ) {
      
      if(is.numeric(bounds)) {
        writeLines(paste("\t L bound : ",bounds[1], sep=""))
        writeLines(paste("\t U bound : ",bounds[2], sep=""))
      }
      
      if( (bounds == "min") | (method == "min") ) {
        writeLines(paste("\t L bound : ",round(min(g),4), sep=""))
        writeLines(paste("\t U bound : ",round(max(g),4), sep=""))
      }
      
    }
    
    writeLines(paste("Mean : ",round(mean(weightsRatio),4), sep=""))
    quantileRatio <- round(quantile(weightsRatio, probs=c(0,0.01,0.1,0.25,0.5,0.75,0.9,0.99,1)),4)
    print(quantileRatio)
    
  }
  
  ## Export in TeX
  if(!is.null(exportDistributionTable)) {
    
    # Linear or raking ratio
    if(is.null(bounds)) {
      
      newNames <- names(quantileRatio)
      newNames <- c(newNames,"Mean")
      
      statsRatio <- c(quantileRatio,mean(quantileRatio))
      names(statsRatio) <- newNames
      
    } else {
      newNames <- names(quantileRatio)
      newNames <- c("L",newNames,"U","Mean")
      
      statsRatio <- c(bounds[1],quantileRatio,bounds[2],mean(quantileRatio))
      names(statsRatio) <- newNames
    }
    
    latexQuantiles <- xtable(as.data.frame(t(statsRatio)))
    
    # Notice that there is one extra column in align(latexQuantiles)
    # since we haven't specified yet to exclide rownames
    if(is.null(bounds)) {
      align(latexQuantiles) <- "|c|ccccccccc||c|"
    } else {
      align(latexQuantiles) <- "|c|c|ccccccccc|c||c|"
    }
    
    
    print(latexQuantiles,  include.rownames = FALSE, include.colnames = TRUE,
           floating = FALSE, file=exportDistributionTable)
  }
  
  
  if(description) {
    writeLines("")
    writeLines("################### Comparison Margins Before/After calibration ###################")
    print(calibrationMarginStats(data=data, marginMatrix=marginMatrix, popTotal=popTotal, colWeights=colWeights, colCalibratedWeights=colCalibratedWeights))
  }
  
  # Plot density of weights ratio
  if(description) {
    if(require("ggplot2")) {
      
      densityPlot = ggplot(data.frame(weightsRatio), aes(x=weightsRatio)) + geom_density(alpha=0.5, fill="#FF6666", size=1.25, adjust=2) + theme_bw()
      print(densityPlot)
      
      if(!is.null(exportDistributionImage)) {
        ggsave(densityPlot, file=exportDistributionImage)
      }
      
    } else {
      warning("Require package ggplot2 to plot weights ratio")
    }
  }
  
  return(g*weights)
}
