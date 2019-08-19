[![Mentioned in Awesome Official Statistics ](https://awesome.re/mentioned-badge.svg)](http://www.awesomeofficialstatistics.org)

# Icarus

Icarus (Icarus Calibrates And Reweights Units in Samples) is an R package providing useful functions for calibration and reweighting estimators in survey sampling. The former name of this package was gaston.

## Cite

To cite Icarus in publications use: Rebecq, Antoine (2017). Icarus: an R package for calibration in survey sampling. R package version 0.3.0.

## Install

You can use the following instruction to install icarus (from CRAN):

```
install.packages("icarus")
```

However, if you wish to install the latest version of icarus, you can use devtools and install directly from this github repo:

```
install.packages("devtools")
library(devtools)
install_github("haroine/icarus")
````

## Short example

In this example, we perform calibration (with the "raking" method) on the test dataset _data_employees_ included in icarus:

```
library(icarus)

N <- 300 ## Population size
## Compute the Horvitz-Thompson estimator (returns 1.666667)
weightedMean(data_employees$movies, data_employees$weight, N)

## Add calibration margins
mar1 <- c("category",3,80,90,60)
mar2 <- c("sex",2,140,90,0)
mar3 <- c("department",2,100,130,0)
mar4 <- c("salary", 0, 470000,0,0)
margins <- rbind(mar1, mar2, mar3, mar4)
## Compute calibration weights
wCal <- calibration(data=data_employees, marginMatrix=margins, colWeights="weight"
                           , method="raking", description=FALSE)
                           
## Value of the calibrated estimator: 2.471917
weightedMean(data_employees$movies, wCal, N)
```
