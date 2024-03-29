% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{dataPop}
\alias{dataPop}
\title{Test population for Icarus.}
\format{
1 column "ident" with unique id for all units.
11 columns with various characteristics of interest for units in the population.
1 column "weight", with sampling weights . Weights equal to zero means that the unit
is not selected in the sample.
1 column "simul_nr" indicates the probability that each unit will respond to the survey.
1 column "responding". For sampled units, indicates whether unit is respondant to survey (1)
or not (0). Variable is also equal to 0 for units not selected in sample
1 column "qTest" containing randomly generated q weights used in unit tests
50000 rows, 1 row per unit in the population.
}
\usage{
dataPop
}
\description{
This data set features a generated population of 50000 units.
11 characteristics of interest for all units in population are featured. 
These characteristics of interest are variously correlated to one another.
A stratified random sampling (with a proportional allocation on variable Y3) 
of fixed size 1000 is selected. Among the 1000 units in the selected sample, only
718 are respondant to the survey. These responding units are selected using a dummy
logit model.
}
\references{
Rebecq, A., & Merly-Alpa, T. Pourquoi minimiser la dispersion des 
poids en sondage. preprint.
}
\author{
Antoine Rebecq
}
