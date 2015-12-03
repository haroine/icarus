
# TODO : better header for file
#############################################
############   Calmar Functions   ###########
#### Functions designed to calibrate  #######
### samples in a familiar INSEE setting #####
#############################################

# Remarque : calmarMatrix = "matrice des marges sans la colonne des noms"

nModalities = function(col)
{
  return(length(unique(col)))
}

calibrationMatrix = function(entryMatrix, popVector=TRUE, isQuantitative=NULL)
{
  if(is.null(isQuantitative)) {
    isQuantitative <- rep(FALSE, ncol(entryMatrix))
  }

  entryMatrix = data.matrix(entryMatrix)

  # Initialization of return matrix
  nRows = nrow(entryMatrix)
  nCols = 0

  N = ncol(entryMatrix)

  # Particular case if entryMatrix has only one row
  if(is.null(N)) {
    N=1
    nRows = length(entryMatrix)
  }

  for(i in 1:N)
  {
    nCols = nCols + nModalities(entryMatrix[,i])
  }

  namesMatrix = names(entryMatrix)
  calibrationMatrix = matrix(0, nRows, 0, byrow=T)
  for(i in 1:N)
  {
    if(!isQuantitative[i]) {
      calibrationMatrix = cbind(calibrationMatrix, colToDummies(entryMatrix[,i], namesMatrix[i]))
    } else {
      calibrationMatrix = cbind(calibrationMatrix, entryMatrix[,i])
    }
  }

  # Add "population" vector
  if(popVector) {
    calibrationMatrix = cbind(calibrationMatrix, rep(1,nrow(calibrationMatrix)))
  }

  return(calibrationMatrix)
}

dummyModalitiesMatrix = function(entryMatrix)
{
  dmatrix = calibrationMatrix(entryMatrix)
  dmatrix[dmatrix!=0] = 1
  return(dmatrix)
}

# TODO : move out of calmarFunctions ?
# (or at least should be "private")
HTtotals = function(dummyModalitiesMatrix, weights)
{
  return(weights%*%dummyModalitiesMatrix)
}

# TODO : ATTENTION, le format change (30/6/14) !!!
# -> maintenant on met le nom de la colonne a gauche
# => "createCalibrationMatrix" est une fonction de compatibilité avec l'ancienne version
createCalibrationMatrix = function(marginMatrix, data, popVector=TRUE)
{
  # Selection des variables de calage dans la table
  # (ainsi que leur caractere qualitatif / quantitatif)
  selectVector = marginMatrix[,1]
  isQuantitative = as.numeric(marginMatrix[,2])

  isQuantitative[isQuantitative != 0] <- 1
  isQuantitative <- 1 - as.numeric(isQuantitative) # is considered as boolean by R

  Xs = data[,selectVector]
  # Mise en forme de la matrice de calage
  matrixCal = calibrationMatrix(Xs, popVector, isQuantitative)

  return(matrixCal)
}

formatMargins = function(calmarMatrix, calibrationMatrix, popTotal=NULL)
{
  # Create empty vector of margins
  cMatrixCopy = calmarMatrix
  if(is.vector(cMatrixCopy)) {
    cMatrixCopy = t(as.matrix(calmarMatrix))
    calmarMatrix = t(as.matrix(calmarMatrix))
  }
  typeMargins = cMatrixCopy[,1]
  typeMargins[typeMargins==0] = 1
  cMargins = rep(0,sum(typeMargins))

  # Fill cMargins
  i=1
  curRow = 1
  while(curRow <= nrow(calmarMatrix))
  {
    if(calmarMatrix[curRow,1] == 0)
    {
      cMargins[i]=calmarMatrix[curRow,2]
      i=i+1
    }
    else
    {
      if(!is.null(popTotal)) {
        popTotalNum <- popTotal
      } else {
        popTotalNum <- 1
      }

      n = calmarMatrix[curRow,1]

      ## If categorial margins are not entered as percentages,
      ## do not multiply by popTotal (except if it is popVector !)
      if( all(calmarMatrix[curRow,2:(n+1)] >= 1) ) {
        popTotalNum <- 1
      }

      for(j in 2:(n+1))
      {
        cMargins[i] = calmarMatrix[curRow,j]*popTotalNum
        i = i+1
      }
    }
    curRow = curRow+1
  }

  # If there is still one column, it is the population one, so we add popTotal to cMargins
  # ... unless specified otherwise
  if(i <= ncol(calibrationMatrix) && !is.null(popTotal))
    cMargins[i] = popTotal

  return(cMargins)
}


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
                       , infinity=1e7, uCostPenalized=1e2, lambda=NULL, gap=NULL
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

    # TODO : Add check marginMatrix -> check if number of modalities in calibration
    # variables matches marginMatrix
    if(!checkNumberMargins(data, marginMatrix)) stop("Error in number of modalities.")

    # TODO : check that sum on all categorial variables is the same

  }

  ## TODO : replace using createFormattedMargins
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
                          , q=rep(1,length(d)), maxIter=maxIter, description=description)
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
          , floating = FALSE, file=exportDistributionTable)
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

#' Gives stats about the calibration process: totals after/before calibration vs. margins
#' Same as first panels output in Calmar/Calmar 2
#' @param data dataframe containing the survey data
#' @param marginMatrix matrix of margins
#' @param popTotal total of population, useful if margins are entered in relative value
#' @param colWeights name of weights column in the dataframe
#' @export
calibrationMarginStats = function(data, marginMatrix, popTotal=NULL, colWeights="POIDS", colCalibratedWeights=NULL, calibThreshold=1.0) {

  displayCalibratedWeights <- TRUE

  if(is.null(colCalibratedWeights)) {
    displayCalibratedWeights <- FALSE
    colCalibratedWeights <- colWeights
  }

  if(displayCalibratedWeights) {
    textAfter <- "After Calibration"
  } else {
    textAfter <- "Current"
  }

  enteredAsPct <- FALSE
  popTotalMarginDisplay <- popTotal
  if(is.null(popTotal)) {
    enteredAsPct <- TRUE
    popTotal <- sum(data[colCalibratedWeights])
    popTotalMarginDisplay <- NA
  }

  toWarn = FALSE
  displayWarningMessage = FALSE

  # Somme des poids (total)
  totalWeights = sum(data.matrix(data[colWeights]))
  totalCalibrated = sum(data[colCalibratedWeights])

  vecTotal = c(totalWeights, totalCalibrated, popTotalMarginDisplay)
  names(vecTotal) = c("Before calibration",textAfter, "Margin")

  vecTotal = round(vecTotal,2)

  marginStatsList = list(vecTotal)

  marginNames = marginMatrix[,1]

  if(is.null(marginMatrix)) {
    names(marginStatsList) = c("Total")
    return(marginStatsList)
  }

  # Other margins
  for(i in 1:nrow(marginMatrix)) {

    toWarn = FALSE
    vecTotal = NULL

    if(as.numeric(marginMatrix[i,2]) == 0) { # If variable is numeric

      sumWeights = data.matrix(data[marginNames[i]])[,1] %*% data.matrix(data[colWeights])[,1]
      sumCalibrated = data.matrix(data[marginNames[i]])[,1] %*% data.matrix(data[colCalibratedWeights])[,1]
      margin = as.numeric(marginMatrix[i,3])

      vecTotal = c(sumWeights, sumCalibrated, margin)
      vecTotal = as.numeric(vecTotal)
      vecTotal = round(vecTotal,2)

      # Check if calibration is exact
      if(is.na(sumCalibrated)) stop(paste("Modality is present in margin tables but not in sample : ",i,";",j))
      if(abs(sumCalibrated - margin) >= calibThreshold) {
        toWarn = TRUE
        displayWarningMessage = TRUE
        #vecTotal = c(vecTotal,"*") # Old convention (same as in Calmar)
        vecTotal = c(vecTotal,round(abs((sumCalibrated - margin)/margin),4))
      }

      if(toWarn == FALSE) {
        names(vecTotal) = c("Before calibration",textAfter,"Margin")
      } else {
        names(vecTotal) = c("Before calibration",textAfter,"Margin", "Warning")
      }
    } else { # If variable has modalities
      modalities = data.matrix(unique(data[marginNames[i]])[,1])
      modalities = sort(modalities)

      # TODO : Assert length(modalities) == marginMatrix[i,2]

      for(j in 1:marginMatrix[i,2]) {

        toWarn = FALSE
        sumWeights = sum(data.matrix(data[data[marginNames[i]] == modalities[j],][colWeights]))
        sumCalibrated = sum(data.matrix(data[data[marginNames[i]] == modalities[j],][colCalibratedWeights]))

        if(!enteredAsPct) {
          margin = as.numeric(marginMatrix[i,2+j])
          tempStatVec = c(sumWeights, sumCalibrated, margin)
        } else {
          margin = as.numeric(marginMatrix[i,2+j])
          tempStatVec = c(sumWeights/totalWeights*100, sumCalibrated/totalCalibrated*100, margin/popTotal*100)
        }

        #tempStatVec = c(sumWeights, sumCalibrated, margin) # TODO : change here level / structure

        # tempStatVec = c(sumWeights/totalWeights*100, sumCalibrated/totalCalibrated*100, margin/popTotal*100)

        tempStatVec = round(tempStatVec,2)


        # Check if calibration is exact
        if(is.na(sumCalibrated)) stop(paste("Modality is present in margin tables but not in sample : ",i,";",j))
        if(abs(sumCalibrated - margin) >= calibThreshold) {
#           toWarn = TRUE
          displayWarningMessage = TRUE
#           tempStatVec = c(tempStatVec, "*")
        }

        vecTotal = rbind(vecTotal, tempStatVec, deparse.level = 0)
      }

      # rownames = marginName_modalities(i)
      rownames(vecTotal) = modalities

      # "Little stars" if not perfectly calibrated
      if(toWarn == FALSE) {
        colnames(vecTotal) = c("Before calibration",textAfter,"Margin")
      } else {
        colnames(vecTotal) = c("Before calibration",textAfter,"Margin", "Warning")
      }



    }

    marginStatsList[[i+1]] = vecTotal
  }


  # Name of statsMargesList
  names(marginStatsList) = c("Total", marginNames)

  if(displayWarningMessage && displayCalibratedWeights)
    writeLines("Careful, calibration may not be exact")

  return(marginStatsList)
}


marginStats <- function(data, marginMatrix, popTotal=NULL, colWeights="POIDS"
                        , colCalibratedWeights=NULL, calibThreshold=1.0) {

  listMarginStats <- calibrationMarginStats(data, marginMatrix, popTotal, colWeights
                                            , colCalibratedWeights, calibThreshold)
  marginStatsDF <- do.call(rbind.data.frame, listMarginStats)

  ## Column difference is re-computed from scratcg
  marginStatsDF <- marginStatsDF[,-c(4)]
  if( is.null(colCalibratedWeights) ) {

    marginStatsDF <- marginStatsDF[,-c(2)] # Do not display calibrated weigths column
    colnames(marginStatsDF) <- c("Before calibration","Margin")
    marginStatsDF$difference <- round(abs(data.matrix(marginStatsDF["Margin"]) - data.matrix(marginStatsDF["Before calibration"]))/data.matrix(marginStatsDF["Margin"])*100,2)

  } else {

    colnames(marginStatsDF) <- c("Before calibration","After calibration","Margin")
    marginStatsDF$difference <- round(abs(data.matrix(marginStatsDF["Margin"]) - data.matrix(marginStatsDF["After calibration"]))/data.matrix(marginStatsDF["Margin"])*100,2)


  }



  return(marginStatsDF)

}


# Does the hot-deck imputation of NAs in calibration variables
# Hot-deck neighbors are selected via margins (of marginMatrix)
## TODO : deprecate and remove
imputCalibrationVars = function(data, marginMatrix) {

  dataUpdated = data
  N = nrow(marginMatrix)

  for(i in 1:N) {

    marginName = marginMatrix[i,1]
    writeLines(paste("Imputation of column : ",marginName))

    if(nrow(data.matrix(data[is.na(data[marginName]),])) > 0) {
      vecParams = marginMatrix[,1]
      vecParams = vecParams[-i]

      dataUpdated[is.na(dataUpdated[marginName]),][marginName] = imputViaNeighbors(data, colToImput=marginName, vecParams, method="first")
    }

  }

  return(dataUpdated)
}


# Check validity of marginMatrix
checkMarginMatrix = function(marginMatrix) {

  checkMatrix = FALSE

  if(is.null(marginMatrix)) return(TRUE) # Case NULL is OK

  # TODO :
  # Check if there are : 1 names column, 1 modalities column and
  # n other columns with n = max(modalities)


  # Check if sum(lines where modalities >=2) = 1.000000


  return(checkMatrix)
}

# Displays number of NAs among margins
missingValuesMargins = function(data, marginMatrix) {

  nVar = nrow(marginMatrix)
  marginNames = marginMatrix[,1]
  returnMatrix = cbind(marginNames, rep(0,nVar))

  for(i in 1:nVar) {
    returnMatrix[i,2] = nrow(data[is.na(data[marginNames[i]]),])
  }

  colnames(returnMatrix) = c("Margin","Missing values")

  return(returnMatrix)
}

# Checks if number of modalities in data matches expected ones according
# to marginMatrix
checkNumberMargins = function(data, marginMatrix) {

  returnBool = TRUE
  marginNames = marginMatrix[,1]

  for(i in 1:length(marginNames)) {

    nModalities = length(table(data.matrix(data[marginNames[i]])))
    expectedModalities = as.numeric(marginMatrix[i,2])
    if(nModalities != expectedModalities && expectedModalities > 0) { ## "0" indicates calibration is made on quantitative total
      writeLines(paste("Error on column ",marginNames[i]," : ",nModalities," modalities in data and ",expectedModalities," expected in margins"))
      return(FALSE)
    }

  }

  return(TRUE)
}


#' Regroups modalities entered in "vecModalities" into single
#' "newModality" in "calibrationMatrix" and adapts "marginMatrix" to the new concept.
#' Typical usage is right before a calibration (and after comptutation of marginMatrix), when
#' you realise calibration output is better when several modalities are reduced to one.
#' (typically very rare modalities, on which calibration constraints are very restrictive).
#' Uses pseudo-"call by reference" via eval.parent because 2 objects are modified :
#' calibrationMatrix and marginMatrix
#' @export
regroupCalibrationModalities <- function(calibrationMatrix, marginMatrix, calibrationVariable, vecModalities, newModality) {

  # First, check if number of modalities match in calibrationMatrix and marginMatrix,
  # otherwise stop
  if(!checkNumberMargins(calibrationMatrix, marginMatrix))
    stop("Number of modalities must match between calibrationMatrix and marginMatrix to regroup calibration modalities.")

  newCalibrationMatrix <- calibrationMatrix
  newMarginMatrix <- marginMatrix

  ## Modification in calibrationMatrix
  newCalibrationMatrix[calibrationVariable] <- regroupUnContiguuousModalities(data.matrix(newCalibrationMatrix[calibrationVariable]), vecModalities, newModality)

  ## Modification in marginMatrix
  # Rank of modified modalities in calibrationMatrix
  calVarModalities <- unique(data.matrix(calibrationMatrix[calibrationVariable]))
  modalitiesRanks <- rank(calVarModalities)
  rankModalitiesDF <- as.data.frame(cbind(calVarModalities,modalitiesRanks))

  vecRanks <- rankModalitiesDF[data.matrix(rankModalitiesDF[calibrationVariable]) %in% vecModalities,]["modalitiesRanks"][,1]
  modifiedLine <- marginMatrix[marginMatrix[,1] == calibrationVariable,]

  # Modify number of modalities in second column of marginMatrix
  # If number of modalities == 1, stop
  oldModalitiesNumber = as.numeric(modifiedLine[2])
  modifiedLine[2] <- as.character(oldModalitiesNumber - length(vecModalities) + 1)

  if( modifiedLine[2] == "1" ) {
    stop("Margin matrix cannot reduce non-numeric margins to numeric margins")
  }

  regroupSum = sum(as.numeric(modifiedLine[vecRanks+2]))

  # Delete regrouped modalities from line
  deletedModalities <- vecRanks
  modifiedLine[deletedModalities+2] <- NA
  modifiedLine <- modifiedLine[!is.na(modifiedLine)]

  newCalVarModalities <- unique(data.matrix(newCalibrationMatrix[calibrationVariable]))
  newModalitiesRanks <- rank(newCalVarModalities)
  rankNewModalitiesDF <- as.data.frame(cbind(newCalVarModalities,newModalitiesRanks))

  insertPosition = rankNewModalitiesDF[rankNewModalitiesDF[,1] == newModality,]$newModalitiesRanks # Rank of new modality in calibrationMatrix modalities

  # Insert regroupSum in modifiedLine
  if(insertPosition+2 > length(modifiedLine)) {
    modifiedLine[insertPosition+2] <- regroupSum
  } else {
    modifiedLine <- c(modifiedLine[1:(insertPosition+2-1)],regroupSum,modifiedLine[(insertPosition+2):length(modifiedLine)])
  }

  # Add 0s to end line
  modifiedLine <- c(modifiedLine, rep("0.0000",ncol(marginMatrix) - length(modifiedLine)))

  # Careful, sum of weights must be equal to 1 even after modalities have been regrouped
  sumMarginLine <- sum(as.numeric(modifiedLine[3:length(modifiedLine)]))

  if( sumMarginLine != 1 ) {
    #print(modifiedLine)
    maxMarginValue <- max(as.numeric(modifiedLine[3:as.numeric(modifiedLine[2])+2]))
    #print(maxMarginValue)
    maxIndex <- which.max(as.numeric(modifiedLine[3:as.numeric(modifiedLine[2])+2]))
    #print(maxIndex)
    modifiedLine[maxIndex+2] <- maxMarginValue + 1 - sumMarginLine
    #print(modifiedLine)
  }

  # Replace in marginMatrix
  newMarginMatrix[marginMatrix[,1] == calibrationVariable,] <- modifiedLine

  # Check if last column of margin matrix is all 0s. If it is, drop last column
  # (means larger line has been reduced). Continue to do so until last colmun is not only 0s.
  while( sum(as.numeric(newMarginMatrix[,ncol(newMarginMatrix)])) == 0 ) {
    newMarginMatrix <- newMarginMatrix[, -ncol(newMarginMatrix)]
  }

  eval.parent(substitute(calibrationMatrix <- newCalibrationMatrix))
  eval.parent(substitute(marginMatrix <- newMarginMatrix))
}


#' Adds a margin to marginMatrix
#'
#' @param marginMatrix The matrix of margins to add the new margin to
#' @param varName Name of variable in calibration matrix corresponding
#' to the new margin
#' @param vecTotals values of margins (Calmar style) for the variable.
#' Note : if length(vecTotals) > 1, then sum(thresholdAdjustToOne) has to be 1.
#' @param adjustToOne if TRUE and sum(vecTotals) is nearly 1, modify values of vecTotals
#' so that sum is 1.
#' @param thresholdAdjustToOne adjust sum(vecTotals) to 1 if difference
#' is under thresholdAdjustToOne
#'
#' @export
addMargin <- function(marginMatrix, varName, vecTotals, adjustToOne=TRUE, thresholdAdjustToOne = 0.01) {

  newMarginMatrix <- marginMatrix

  # Length of vecTotals :
  if( length(vecTotals) == 1 ) {
    nModality <- 0
  } else {
    if( length(vecTotals) > 1 ) {
      nModality <- length(vecTotals)
    } else {
      stop("vecTotals must be non NULL vector")
    }
  }

  # TODO : adjust vecTotals to 1
  if( nModality > 1 && sum(vecTotals) != 1  ) {

    if(adjustToOne && abs(sum(vecTotals) - 1) < thresholdAdjustToOne) {
      # TODO : adjust highest value
      maxMarginValue <- max(as.numeric(vecTotals))
      maxIndex <- which.max(as.numeric(vecTotals))
      vecTotals[maxIndex] <- maxMarginValue + 1 - sum(vecTotals)
    } else {
      stop("sum(vecTotals) must be equal to 1.")
    }

  }

  newMarginLine <- c(varName, nModality, vecTotals)

  # newMarginLine must have right format before it is added to
  # newMarginMatrix
  if(length(newMarginLine) < ncol(newMarginMatrix)) {
    # Add missing zeroes :
    missingZeroes <- rep(0, ncol(newMarginMatrix) - length(newMarginLine))
    newMarginLine <- c(newMarginLine, missingZeroes)
  }

  if(length(newMarginLine) > ncol(newMarginMatrix)) {
    # Add columns of 0s to newMarginMatrix
    missingZeroes <- matrix(nrow = nrow(newMarginMatrix), ncol = (length(newMarginLine) - ncol(newMarginMatrix)), 0)
    newMarginMatrix <- cbind(newMarginMatrix, missingZeroes)
  }

  # Append to newMarginMatrix :
  newMarginMatrix <- rbind(newMarginMatrix, newMarginLine, deparse.level = 0)

  return(newMarginMatrix)
}

## Modifies margin
modifyMargin <- function(marginMatrix, varName, vecTotals, adjustToOne=TRUE, thresholdAdjustToOne = 0.01) {

  # Delete selected margin
  indexSelectedMargin <- NULL
  i <- 1
  while(i <= nrow(marginMatrix)) {
    if(marginMatrix[i,1] == varName) {
      indexSelectedMargin <- i
    }
    i <- i+1
  }

  newMarginMatrix <- marginMatrix[-indexSelectedMargin,]
  if(is.null(ncol(newMarginMatrix))) {
    newMarginMatrix <- t(as.matrix(newMarginMatrix))
  }

  # Add selected margin
  newMarginMatrix <- addMargin(newMarginMatrix, varName, vecTotals, adjustToOne, thresholdAdjustToOne)

  return(newMarginMatrix)
}

## Private function that creates margins to the right format
createFormattedMargins <- function(data, marginMatrix, popTotal=NULL) {

  if(is.null(marginMatrix)) {

    if(is.null(popTotal)){
      stop("No margin or population total specified for dataMen.")
    }

    writeLines("Calibration only made on population totals for dataMen")
    matrixCal = rep(1,nrow(data))
    formattedMargins = c(popTotal)

  } else {

    # Creation of the elements
    calmarMatrix = marginMatrix[,2:ncol(marginMatrix)]
    # Transform calmarMatrix to numeric matrix to avoid problems in formatMargins
    if(!is.vector(calmarMatrix)) {
      calmarMatrix = matrix(as.numeric(calmarMatrix), nrow=nrow(calmarMatrix), ncol=ncol(calmarMatrix), byrow=F)
    } else {
      calmarMatrix = as.numeric(calmarMatrix)
    }

    popVector <- TRUE
    if(is.null(popTotal)) {
      popVector <- FALSE
    }

    matrixCal = createCalibrationMatrix(marginMatrix,data, popVector)

    formattedMargins = formatMargins(calmarMatrix, matrixCal, popTotal)

  }

  return(list(formattedMargins, matrixCal))

}

## TODO : documentation about integrated calibration
integratedCalibration <- function(dataMen, marginMatrixMen, popMen = NULL,
                                    dataInd, marginMatrixInd, popInd = NULL,
                                    identMen="IDENT_LOG", identInd="IDENT_IND",
                                    colWeights = "POIDS", colCalibratedWeights="POIDS_CALES", method="linear",
                                    maxIter=2500, description=FALSE, bounds, scale=TRUE, check=TRUE) {

  ## TODO : check that idents given are present in tables

  ## Merge right, by identMen
  dataSimultaneous <- merge(dataMen, dataInd, by=identMen, all.y=TRUE)

  #### Count number of ind margin variables by men unit
  # For continuuous variables
  quantiIndMargins <- marginMatrixInd[marginMatrixInd[,2]=="0",]
  if(is.null(ncol(quantiIndMargins))) {
    quantiIndMargins <- quantiIndMargins[1]
  } else {
    quantiIndMargins <- quantiIndMargins[,1]
  }
  print(quantiIndMargins) # debugging
  # For discrete variables
  qualiIndMargins <- marginMatrixInd[marginMatrixInd[,2]!="0",]
  if(is.null(ncol(qualiIndMargins))) {
    qualiIndMargins <- qualiIndMargins[1]
  } else {
    qualiIndMargins <- qualiIndMargins[,1]
  }
  print(qualiIndMargins) # debugging

  ## TODO : quanti variables to dummies
  ## TODO : aggregate individual variables by identMen

  ### Format margins in one margin table
  # Margins for dataMen : TODO -> remove
#   formattedMarginsMen <- createFormattedMargins(dataMen, marginMatrixMen, popMen)[[1]]
#   matrixCalMen <- createFormattedMargins(dataMen, marginMatrixMen, popMen)[[2]]
#   # Margins for dataInd
#   formattedMarginsInd <- createFormattedMargins(dataInd, marginMatrixInd, popInd)[[1]]
#   matrixCalInd <- createFormattedMargins(dataInd, marginMatrixInd, popInd)[[2]]
#
  ## Add individual margins to marginMatrixMen as quantitative variables
  ## TODO


  #### Calibration

  return(dataSimultaneous)
}

##
