# copyright (C) 2014-2016 A.Rebecq

# TODO : add xtable parameters to marginsToTeX parameters
marginsToTeX <- function(marginMatrix, names=NULL, pct=FALSE, popTotal=NULL,
                         scaleboxTeX=NULL, file=NULL,
                         label=NULL, caption=NULL) {
  
  if (!requireNamespace("xtable", quietly = TRUE)) {
    stop("Package xtable needed for export of margins in LateX to work. Please install it.",
         call. = FALSE)
  }
  
  if(!is.matrix(marginMatrix)) {
    stop("marginsToTeX input type has to be matrix.")
  }
  
  if(!is.null(names)) {
    if(length(names) != nrow(marginMatrix)) {
      stop("Name length must equal number of rows in marginMatrix")
    }
    
    marginMatrix[,1] <- names
  }
  
  if(pct) {
    numericPart <- marginMatrix[,3:(ncol(marginMatrix))]
    numericPart <- as.numeric(numericPart)
    numericPart <- 100*numericPart
    numericPart -> marginMatrix[,3:(ncol(marginMatrix))]
  }
  
  # Write zeros as NA
  marginMatrix[as.numeric(marginMatrix) == 0] <- NA

  marginDF <- as.data.frame(marginMatrix)
  
  # Heuristic rule for scalebox
  if(is.null(scaleboxTeX)) {
    if(ncol(marginDF) >= 10) {
      scaleboxTeX <- 1.4 - ncol(marginDF) / 20
    }
    
    if(ncol(marginDF) >= 28) {
      stop("Automatic scaleboxing not configured for more than 28 margins.")
    } 
  }
  
  captionTeX <- caption
  if(!is.null(popTotal)) {
    captionTeX <- paste(caption, " -- total population : ", round(popTotal,0),sep="")
  }
  
  print(xtable::xtable(marginDF, caption=captionTeX, label=label), include.rownames = FALSE, include.colnames = FALSE,
               floating = TRUE, scalebox=scaleboxTeX, file=file
        )

}

#' Create empty margin matrix
#' @description 
#' Use this to create an empty margin matrix (which facilitates
#' the use of magrittr syntax to enter margins)
#' 
#' @examples 
#' library(magrittr)
#' N <- 230 ## population total
#' ## Horvitz Thompson estimator of the mean: 2.174
#' weightedMean(data_employees$movies, data_employees$weight, N)
#' ## Enter calibration margins:
#' margins <- newMarginMatrix() %>%
#'   addMargin("category", c(0.35, 0.40, 0.25)) %>%
#'   addMargin("sex", c(0.6, 0.4)) %>%
#'   addMargin("department", c(0.45, 0.55)) %>%
#'   addMargin("salary", 470000)
#' ## Compute calibrated weights with raking ratio method
#' wCal <- calibration(data=data_employees, marginMatrix=margins, colWeights="weight"
#'                     , method="raking", pct = TRUE, description=FALSE
#'                     , popTotal = N)
#' ## Calibrated estimate: 2.471917
#' weightedMean(data_employees$movies, wCal, N)
#' @export
newMarginMatrix <- function() {
  return(matrix(, nrow = 0, ncol = 1))
}