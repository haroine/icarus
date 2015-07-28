
#' Create Non-response Adjusted Weights using the HRG method
#'
#' @param data The dataframe containing the survey results. Must include
#' every variable needed by modelGRH as well as a column containing the sampling
#' weights for each responding unit
#' @param modelGRH The logit model predicting the non-response
#' @param colPoids The column name of the column containing the weights
#' in the dataframe
#' @param colNew Specifies the name of the columns added to the dataframe
#' @param nGRH Number of groups created. If NULL, the number of groups is computed
#' using the Beaumont-Haziza method
#' @param If TRUE, display a few stats and graphs about the NRA weights created
#'
#' @return Input dataframe augmented with three columns : numGRH (id of HRG the unit is assigned to),
#' rapportPoidsCNR (the non-response adjustment factor) and
#' POIDS_CNR (the final non-response adjusted weights)
#'
#' @export
ajouterPoidsGRH = function(data, modelGRH, colPoids="POIDS", colNew=c("pHat","numGRH","rapportPoidsCNR","POIDS_CNR"),
                    nGRH=NULL, stats=FALSE) {

  # On recupere les noms des colonnes :
  if(length(colNew) != 4) stop("colNew must have 4 parameters.")
  colpHat = colNew[1]
  numGRH = colNew[2]
  rapportPoidsCNR = colNew[3]
  POIDS_CNR = colNew[4]

  # Ajouter les pHat a la matrice des donnees
  pHat = as.matrix(modelGRH$fitted.values)
  colnames(pHat) = c(colpHat)
  data = merge(data, pHat, by=c("row.names"), all.x=TRUE)


  # Imputer des pHat aux pHat NA (s'il y a lieu) et trier la matrice des donnees par pHat
  if(nrow(data[is.na(data[colpHat]),]) > 0) data[is.na(data[colpHat]),][colpHat] = imputpHat(data, modelGRH=modelGRH)
  data = data[order(data[colpHat]),]

  # Si aucun nombre de GRH est indique, on le calcule par la methode de Haziza-Beaumont
  if(!is.numeric(nGRH)) {
    writeLines("No nGRH entered, computing with Haziza-Beaumont method")
    nGRH = ngrhHazizaBeaumont(100, data=data, seuil=0.99)
  }
  data[numGRH] = attribuerGRH(data, method="quantiles", nGRH)
  data[rapportPoidsCNR] = rapportPoidsCNR(data, colPoids=poidsInitial)
  data[POIDS_CNR] = poidsCNR(data, colPoids=poidsInitial)


  # Stats
  if(stats) {
    # Sur les GRH
    statsGRH(data, colRapportPoids="rapportPoidsCNR", colPoids=colPoids)

    # TODO : ajouter selection judicieuse de stats ??
    # (notamment sur les taux de collecte ?)
  }

  return(data)

}


# TODO : documenter la fonction
# TODO : very slow !!... -> improve (C ?)
# Switch from String methods to function methods
## TODO : remove from package
imputViaNeighbors = function(data, colToImput, vecParams, method="mean") {

  # Select NAs from colToImput
  colimputNA = data.matrix(data[is.na(data[colToImput]),])

  returnVector = rep(0,nrow(colimputNA))

  for(i in 1:nrow(colimputNA)) {

    # Fill vecVal with value of row i
    vecVal = rep(0,length(vecParams))
    for(j in 1:length(vecParams)) {

      ident = as.numeric(colimputNA[i,1]) # row.names is always first column
      vecValue = data[data[,1]==ident,][vecParams[j]][1,1] # row.names is always first column

      if(!is.na(vecValue))
        vecVal[j] = vecValue
      else
      {
        # S'il reste au moins deux colonnes dans vecParams, virer la colonne defectueuses
        # (et afficher un warning). S'il reste une seule colonne -> stop.
        if(length(vecParams) >= 2) {

          newVecParams = vecParams[!vecParams %in% c(vecParams[j])]
          writeLines(paste("NA found in colmun. Starting over without column",vecParams[j]))
          return(imputViaNeighbors(data, colToImput, vecParams=newVecParams, method))
        } else {
          stop(paste("Note : NA found in value number ",j," of row number ",i," and no replacement column available"))
        }
      }

    }

    # TODO : penser à mettre un warning si aucun "plus proche voisin" n'est trouvé
    dataVoisins = matrixVoisins(data, vecParams, vecVal)

    if(nrow(dataVoisins) == 0)
      stop("WARNING : no neighbors for selected condition")

    if(nrow(dataVoisins[!is.na(dataVoisins[colToImput]),]) == 0) {
      k = 1 # Start over without first column (arbitrary)
      newVecParams = vecParams[!vecParams %in% c(vecParams[k])]
      writeLines(paste("No non-NA neighbor. Starting over without column",vecParams[k]))
      return(imputViaNeighbors(data, colToImput, vecParams=newVecParams, method))
    }

    switch(method,
      mean={
        imputValue = base::mean(data.matrix(dataVoisins[colToImput]), na.rm=TRUE) # Imput with mean of the "neighbors"
      },
      first={
        imputValue = firstNotNA(data.matrix(dataVoisins[colToImput])) # Imput with value of the first "neighbor"
      },
      median={
        imputValue = base::median(data.matrix(dataVoisins[colToImput]), na.rm=TRUE) # Imput with median of the "neighbors"
      },
      {
        imputValue = base::mean(data.matrix(dataVoisins[colToImput]), na.rm=TRUE) # Imput with mean of the "neighbors"
        print('Default method : mean')
      }
    )

    returnVector[i] = imputValue

  }

  return(returnVector)
}

## TODO : remove from package
firstNotNA = function(vec, i=1)
{
  if(!is.na(vec[i])) {
    return(vec[i])
  } else {
    if(length(vec) >= i)
      firstNotNA(vec, i+1)
    else
      stop("No non-NA value")
  }
}

# TODO : documenter
# Fonction d'imputation des pHat pour les NA
# TODO : par défaut, vecParams = vec du modèle GRH s'ils existe :
# vecParams=names(modelGRH$coefficients)[2:length(modelGRH$coefficients)]
## TODO : imput by only using the mean of the pHat column
imputpHat = function(data, vecParams=NULL, modelGRH=NULL) {
  if(is.null(vecParams)) {
    if(!is.null(modelGRH)) vecParams = names(modelGRH$coefficients)[2:length(modelGRH$coefficients)]
    else stop("Need to enter vecParams or modelGRH")
  }

  return(imputViaNeighbors(data, colToImput="pHat", vecParams))
}


# TODO :  documenter la fonction
# Exemple : test2 = matrixVoisins(enlNPdC, c(enlNPdC$LOGEMENTCO, enlNPdC$LOGEMENT_PETIT), c(1,1))
# TODO : enlever la liste et remplacer par un vecteur de noms de colonnes ? (si possible)
matrixVoisins = function(data, vecParams, vecVal) {
  dataReturn = data

  for(i in 1:length(vecParams))
  {
    dataReturn = dataReturn[dataReturn[vecParams[i]] == vecVal[i],]
  }

  return(dataReturn)
}


# TODO : autres méthodes que quantiles
# Suppose que data est ordonnée selon la colonne des pHat
attribuerGRH = function(data, method="quantiles", nGRH, colRepondant="repondant") {

  numGRH = NULL
  if(method=="quantiles")
  {
    reste = nrow(data)%%nGRH
    sizeGRH = rep(nrow(data)%/%nGRH,nGRH)

    if(reste >= 1)
      for(i in 1:reste) sizeGRH[i] = sizeGRH[i]+1

    for(j in 1:nGRH) numGRH = c(numGRH,rep(j,sizeGRH[j]))
  }

  return(numGRH)
}

rapportPoidsCNR = function(data, colNumGRH="numGRH", colPoids="POIDS", colRepondant="repondant") {

  rapportPoids = NULL

  numGRH = data.matrix(data[colNumGRH])
  nGRH = length(unique(numGRH))
  poidsInit = data.matrix(data[colPoids])

  for(i in 1:nGRH) {
    repondantsInGRH = data[data[colRepondant]==1 & data[colNumGRH]==i,]

    if(nrow(repondantsInGRH) > 0) { # Handle case when there is no unit with colRepondant==1 in GRH
      rowInGRH = data[data[colNumGRH]==i,]

      sommePoidsRepondants = sum(repondantsInGRH[colPoids])
      sommePoids = sum(rowInGRH[colPoids])

      rapportPoids = c(rapportPoids, rep(sommePoids/sommePoidsRepondants,nrow(rowInGRH)))
    } else {
      writeLines(paste("No unit in GRH number : ",i, sep=""))
      rapportPoids = c(rapportPoids, rep(NA,nrow(rowInGRH)))
    }
  }

  return(rapportPoids)
}

poidsCNR = function(data, colNumGRH="numGRH", colPoids="POIDS", colRepondant="repondant") {

  rapport = rapportPoidsCNR(data, colNumGRH, colPoids, colRepondant)

  poidsFinal = data.matrix(data[colPoids])*rapport

  return(poidsFinal)
}

statsGRH = function(data, colRepondant="repondant", colNumGRH="numGRH", colPoids="POIDS", colPoidsCNR="POIDS_CNR", colRapportPoids="",
                    exportPath=NULL, suffixFile="") {

  if(colRapportPoids!="")
    rapport = data.matrix(data[colRapportPoids])
  else
    rapport = (data.matrix(data[colPoidsCNR])) / (data.matrix(data[colPoids]))

  nGRH = length(unique(data.matrix(data[colNumGRH])))

  individusParGRH = rep(0, nGRH)
  repondantsParGRH = rep(0,nGRH)

  for(i in 1:nGRH) {
    repondantsInGRH = data[data[colRepondant]==1 & data[colNumGRH]==i,]
    rowInGRH = data[data[colNumGRH]==i,]
    individusParGRH[i] = nrow(rowInGRH)
    repondantsParGRH[i] = nrow(repondantsInGRH)
  }

  ## Plot histograms

  if(require("ggplot2")) {

    plotRapport <- qplot(as.data.frame(rapport)$rapport, geom="histogram", main="Ratio NRA weights / initial weights", xlab="Ratio", ylab = "Frequency")
    plotIndividus <- qplot(as.data.frame(individusParGRH)$individusParGRH, geom="histogram", main="Units per group", xlab="Number of units", ylab = "Frequency")
    plotRepondants <- qplot(as.data.frame(repondantsParGRH)$repondantsParGRH, geom="histogram", main="Answering units per group", xlab="Number of units", ylab = "Frequency")

    print(plotRapport)
    print(plotIndividus)
    print(plotRepondants)

    # Save images
    if(!is.null(exportPath)) {
      ggsave(plotRapport, filename=paste(exportPath,"ratioNRA_",suffixFile,".pdf",sep=""))
      ggsave(plotIndividus, filename=paste(exportPath,"unitsNRA_",suffixFile,".pdf",sep=""))
      ggsave(plotRepondants, filename=paste(exportPath,"answ_units_NRA_",suffixFile,".pdf",sep=""))
    }

  } else {

    hist(rapport)
    #print("Nombre d'individus par GRH : ")
    hist(individusParGRH)

    #print("Nombre de répondants par GRH : ")
    hist(repondantsParGRH)

  }

}


# Methode de Haziza et Beaumont pour déterminer le nombre de GRH
# TODO : utiliser dichotomie plutôt que for
#' Computes the optimal number of HRG according to the Beaumont-Haziza method
#' @export
ngrhHazizaBeaumont = function(nGRHtests, data, seuil=0.99, colpHat="pHat") {

  rsquaredvec = rep(0, nGRHtests)

  for(nGRH in 1:nGRHtests) {

    numGRH = attribuerGRH(data, method="quantiles", nGRH)

    # Méthode de Beaumont et Haziza :
    dummiesNumGRH = colToDummies(numGRH, "numGRH")

    #regressionLin = lm(data$pHat ~ dummiesNumGRH)
    regressionLin = lm(data.matrix(data[colpHat]) ~ dummiesNumGRH)
    rsquared = summary(regressionLin)$r.squared
    #print(paste("Nombre de GRH : ",nGRH,"; R² = ",rsquared))
    rsquaredvec[nGRH] = rsquared
  }

  nGRHoptimal = nGRHtests-length(rsquaredvec[rsquaredvec>=seuil])+1
  # TODO : warning si le seuil n'est pas atteint

  return(nGRHoptimal)
}







# TODO : documenter (donner les exemples)
# Exemples :
# statsTauxCollecte(enlNPdC, nameCol="OCC_STATUT_OCC", nameDummy="STOCD")
# statsTauxCollecte(enlNPdC, nameCol="TUU", selection=c("PC","GC","RURAL"))
#' Outputs stats about the collection rate
#' @export
statsTauxCollecte = function(data, nameCol, nameDummy=nameCol, colRepondant = "repondant", selection=c("auto"), sepDummies="_") {

  if(selection[1]=="auto")
    modalities = unique(data.matrix(data[nameCol]))
  else
    modalities = selection

  nModalities = length(modalities)
  occurences = rep(0,nModalities)
  tauxCollecte = rep(0,nModalities)

  for(i in 1:nModalities)
  {
    dummyName = paste(nameDummy,sepDummies,modalities[i], sep="")
    occurences[i] = nrow(data[data[dummyName]==1,])
    tauxCollecte[i] = nrow(data[data[dummyName]==1 & data[colRepondant]==1,])/occurences[i]
  }

  statsMatrix = cbind(modalities, occurences, tauxCollecte)
  statsMatrix = statsMatrix[order(statsMatrix[,1]),]

  return(statsMatrix)
}






# TODO : documenter
# TODO : par défaut, vecMarges devrait correspondre aux colonnes de la table de marges si elle existe
# TODO : Si la table des marges existe, ajouter la colonne marge à côté
#' Outputs stats about the non-response adjustment process for some variables
#' @param vecMarges names of variables for which stats are computed
#' @return List of stats for each vector
#' @export
statsMarges = function(data, vecMarges, colPoids = "POIDS", colPoidsCNR="POIDS_CNR", colRepondant="repondant", sepDummies="_", statsEchantillon=FALSE) {

  # Somme des poids (total)
  totalEchantillon = sum(data.matrix(data[colPoids]))
  totalPoids = sum(data.matrix(data[data[colRepondant]==1,][colPoids]))
  totalCNR = sum(data.matrix(data[data[colRepondant]==1,][colPoidsCNR]))

  vecTotal = c(totalPoids, totalCNR)
  names(vecTotal) = c("Avant CNR","Après CNR")

  if(statsEchantillon) {
    namesSave = names(vecTotal)
    vecTotal = c(totalEchantillon, vecTotal)
    names(vecTotal) = c("Echantillon",namesSave)
  }

  vecTotal = round(vecTotal,0)

  statsMargesList = list(vecTotal)

  for(i in 1:length(vecMarges)) {

    colName = vecMarges[i]
    dummyName = paste(colName, sepDummies,sep="")
    # Check if there are dummies associated with margin
    dummies = grepl(dummyName,colnames(data))
    dummies = colnames(data)[dummies==TRUE]
    # Order alphabetically dummies
    dummies = sort(dummies)
    # TODO : careful, "dummyNames_other" might not be last column

    if(length(dummies)==0)
    {
      statMarge = c(sum(data.matrix(data[data[colName]!=0 & data[colRepondant]==1,][colPoids])),
                    sum(data.matrix(data[data[colName]!=0 & data[colRepondant]==1,][colPoidsCNR])))

      statMarge = round(statMarge,0)

      names(statMarge) = c("Avant CNR","Après CNR")

      if(statsEchantillon) {
        namesSave = names(statMarge)
        statMarge = c(sum(data.matrix(data[data[colName]!=0,][colPoids])),
                      statMarge)
        names(statMarge) = c("Echantillon",namesSave)
      }
    }
    else
    {

      statMargeTotal0 = NULL
      statMargePourcentage0 = NULL
      statMargeTotal1 = NULL
      statMargePourcentage1 = NULL
      statMargeTotal2 = NULL
      statMargePourcentage2 = NULL

      for(j in 1:length(dummies)) {

        sommePoids0 = sum(data.matrix(data[data[dummies[j]]==1,][colPoids]))
        sommePoids1 = sum(data.matrix(data[data[dummies[j]]==1 & data[colRepondant]==1,][colPoids]))
        sommePoids2 = sum(data.matrix(data[data[dummies[j]]==1 & data[colRepondant]==1,][colPoidsCNR]))

        statMargeTotal0 = c(statMargeTotal0,
                            sommePoids0
        )
        statMargePourcentage0 = c(statMargePourcentage0,
                                  sommePoids0/totalEchantillon*100
        )

        statMargeTotal1 = c(statMargeTotal1,
                            sommePoids1
        )
        statMargePourcentage1 = c(statMargePourcentage1,
                                  sommePoids1/totalPoids*100
        )

        statMargeTotal2 = c(statMargeTotal2,
                            sommePoids2
        )
        statMargePourcentage2 = c(statMargePourcentage2,
                                  sommePoids2/totalCNR*100
        )

      }

      statMarge = rbind(statMargeTotal1,statMargePourcentage1,statMargeTotal2,statMargePourcentage2)

      colnames(statMarge) = dummies
      rownames(statMarge) = c("Total avant CNR", "Pourcentage avant CNR","Total après CNR", "Pourcentage après CNR")

      statMarge[1,] = round(statMarge[1,],0)
      statMarge[2,] = round(statMarge[2,],2)
      statMarge[3,] = round(statMarge[3,],0)
      statMarge[4,] = round(statMarge[4,],2)

      if(statsEchantillon) {
        rowNamesSave = rownames(statMarge)
        statMarge = rbind(statMargeTotal0,statMargePourcentage0,statMarge)
        rownames(statMarge) = c("Total échantillon", "Pourcentage échantillon", rowNamesSave)


        statMarge[1,] = round(statMarge[1,],0)
        statMarge[2,] = round(statMarge[2,],2)
      }
    }

    statsMargesList[[i+1]] = statMarge
  }

  # name of statsMargesList
  names(statsMargesList) = c("Total", vecMarges)

  return(statsMargesList)
}
