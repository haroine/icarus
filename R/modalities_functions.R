# TODO : Mettre ces fonctions dans un fichier source séparé

# TODO : recodeModality function

# TODO : regroupModalities function : pour faire des tranches

#' Changes a column containing multiple values
#' to a matrix of columns containing the dummies
#' corresponding to each value.
#'
#' @param col input column
#' @param nameCol  name that will be used as a prefix for
#' dummies column name in the output matrix
#' @param modalities if a vector is entered, dummies produced
#' will only be the ones corresponding to the values in the
#' "modalities" input column + another one containing all
#' the other modalities.
#' @param keepValue Logical. If TRUE, puts not "1"s in the dummies output
#' columns but the real values in the "col" column (except if values are
#' non-numeric)
#' @return Matrix containing the dummy columns
#' @export
colToDummies = function(col, nameCol, modalities=NULL, keepValue = FALSE)
{
  if(is.null(modalities))
  {
    modalities = unique(col)
    fillWithRest = FALSE
  }
  else
  {
    modalities = unique(modalities) # In case there are repeated values in input modalities vector
    fillWithRest = TRUE
  }

  nModalities = length(modalities)

  # Create a matrix with as many vectors as there are modalities
  if(fillWithRest)
    dummyMatrix = matrix(0, length(col), nModalities+1, byrow=T)
  else
    dummyMatrix = matrix(0, length(col), nModalities, byrow=T)
  names = rep(NULL,nModalities)

  # Fill with dummies
  N = nModalities
  for(i in 1:N)
  {
    if(keepValue && is.numeric(modality))
      modality = modalities[i]
    else
      modality = 1

    dummyMatrix[,i][col==modalities[i]] = modality
    names[i] = paste(nameCol,modalities[i],sep="_")
  }

  # Names are "originalname"_modality
  if(fillWithRest) names[N+1] = "ZZZZZZZZZ" # Ugly hack useful to sort by colnames and leaving "modality_other" at the end

  colnames(dummyMatrix) = names

  # order columns by alphanumeric order
  dummyMatrix = dummyMatrix[,order(colnames(dummyMatrix))]

  if(fillWithRest)
  {
    dummyMatrix[,nModalities+1][!(col %in% modalities)] = 1
    names[nModalities+1] = paste(nameCol,"other",sep="_")
    colnames(dummyMatrix)[nModalities+1] =names[nModalities+1]
  }

  return(dummyMatrix)
}


###############################################
###############################################

# "modalities" = replacement Modalities
# If "modalities" is not specified, replacement modalities will be 1 -> length(column)
# TODO : add option for "modalities" to keep first modality of each group as replacement modality
# Example : regroupModalities(enlDataNeufs$AGEPR, regroupMatrix)
#' Regroup
#' @param column Column vector which values are going to be replaced
#' @param regroupMatrix Bounds of the values to regroup under the same modality
#' @param modalities Specify the values of the modalities to use. Must match number of rows
#' of regroupMatrix
#' If not specified, replacement modalities will be 1:length(column)
#' @return Column vector with regrouped modalities
#' @export
#' ## TODO: example
regroupModalities = function(column, regroupMatrix, modalities=NULL) {

  regroupedColumn = rep(NA, length(column))

  if(is.null(modalities))
    modalities = c(1:nrow(regroupMatrix))

  for(i in 1:nrow(regroupMatrix)) {
    regroupedColumn[column >= regroupMatrix[i,1] & column <= regroupMatrix[i,2]] = modalities[i]
  }

  return(regroupedColumn)

}

# As regroupModalities is designed to regroup continuuous modalities
# (typically to transform numeric variable into discrete ones),
# regroupContiguuousModalities can regroup different modalities into just one
regroupUnContiguuousModalities <- function(column, vecModalities, newModality) {

  regroupedColumn <- column

  regroupedColumn[regroupedColumn %in% vecModalities] = newModality

  return(regroupedColumn)

}
