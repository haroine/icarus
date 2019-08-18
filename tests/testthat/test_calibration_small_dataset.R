# copyright (C) 2014-2016 A.Rebecq
library(testthat)

context("Test calibration functions on small dataset")

test_that("Calibration functions check out with Calmar", {

  ## Calibration margins
  mar1 <- c("category",3,80,90,60)
  mar2 <- c("sex",2,140,90,0)
  mar3 <- c("department",2,100,130,0)
  mar4 <- c("salary", 0, 470000,0,0)
  margins <- rbind(mar1, mar2, mar3, mar4)


  wCalesLin <- calibration(data=data_employees, marginMatrix=margins, colWeights="weight"
                           , method="linear", description=FALSE)

  wCalesRaking <- calibration(data=data_employees, marginMatrix=margins, colWeights="weight"
                           , method="raking", description=FALSE)

  wCalesLogit1 <- calibration(data=data_employees, marginMatrix=margins, colWeights="weight"
                              , method="logit", description=FALSE, bounds=c(0.4,2.2))

  ## TODO : test M=3 in Calmar

  expect_equal(wCalesLin, calWeights_movies$wLinear, tolerance=1e-5)
  expect_equal(wCalesRaking, calWeights_movies$wRaking, tolerance=1e-5)
  expect_equal(wCalesLogit1, calWeights_movies$wLogit, tolerance=1e-5)

  ## Test estimators value
  expect_equal(HTmean(data_employees$movies, wCalesLin), 2.93, tolerance=1e-2)
  expect_equal(HTmean(data_employees$movies, wCalesRaking), 3.22, tolerance=1e-2)
  expect_equal(HTmean(data_employees$movies, wCalesLogit1), 3.14, tolerance=1e-2)
  
  ## Test errors in convergence
  expect_error(
    calibration(data=data_employees, marginMatrix=margins, colWeights="weight"
              , method="linear", description=FALSE, popTotal = 300),
    "no convergence", ignore.case=T)
  
  expect_error(
    calibration(data=data_employees, marginMatrix=margins, colWeights="weight"
                , method="raking", description=FALSE, popTotal = 300),
    "no convergence", ignore.case=T)
  
  expect_error(
    calibration(data=data_employees, marginMatrix=margins, colWeights="weight"
                , method="logit", bounds=c(0.5,1.5), description=FALSE, popTotal = 300),
    "no convergence", ignore.case=T)
  
  expect_error(
    calibration(data=data_employees, marginMatrix=margins, colWeights="weight"
                , method="logit", bounds=c(0.6,2.0), description=FALSE),
    "no convergence", ignore.case=T)
  
  ## Errors in weights column
  data_employees$weight_bug <- data_employees$weight
  data_employees$weight_bug[c(1,10)] <- c(0,NA)
  
  expect_error(
    calibration(data=data_employees, marginMatrix=margins, colWeights="weight_bug"
              , method="linear", description=TRUE),
    "Some weights are NA", ignore.case=T)
  
  data_employees$weight_bug[c(1,10)] <- c(0,0)
  
  expect_error(
    calibration(data=data_employees, marginMatrix=margins, colWeights="weight_bug"
                , method="linear", description=TRUE),
    "Some weights are negative or zero", ignore.case=T)
  
  
})

test_that("Penalized calibration checks out", {

  mar1 <- c("categ",3,80,90,60)
  mar2 <- c("sexe",2,140,90,0)
  mar3 <- c("service",2,100,130,0)
  mar4 <- c("salaire", 0, 470000,0,0)
  margins <- rbind(mar1, mar2, mar3, mar4)

  costsInfty <- c(Inf, Inf, Inf, Inf)

  expect_warning(
    wCalesLin <- calibration(data=data_employees, marginMatrix=margins, colWeights="weight"
                           , description=FALSE, costs=costsInfty),
    "all costs are infinite", ignore.case=T)

  expect_equal(wCalesLin, calWeights_movies$wLinear, tolerance=1e-4)

  ## Test bad specification of infinite costs
  costsInfty2 <- c(-3, Inf, -1, -5000)
  expect_warning(
    wCalesLin2 <- calibration(data=data_employees, marginMatrix=margins, colWeights="weight"
                              , description=FALSE, costs=costsInfty2),
    "all costs are infinite", ignore.case=T)

  expect_equal(wCalesLin2, wCalesLin, tolerance=1e-7)

  ## TODO: test if too wide gap returns linear calibration
#   wPenalized1 <- calibration(data=data_employees, marginMatrix=margins, colWeights="weight"
#                              , description=TRUE, costs=c(100,10,1,0.1,Inf)
#                              , gap=2, popTotal=230)
#
#   print(wPenalized1)
#   expect_equal(wPenalized1, calWeights_movies$wLinear, tolerance=1e-7)

})


test_that("Calibration on minimum bounds checks out", {

  mar1 <- c("category",3,80,90,60)
  mar2 <- c("sex",2,140,90,0)
  mar3 <- c("department",2,100,130,0)
  mar4 <- c("salary", 0, 470000,0,0)
  margins <- rbind(mar1, mar2, mar3, mar4)

  wCalesMin <- calibration(data=data_employees, marginMatrix=margins, colWeights="weight"
                           , description=FALSE, bounds="min", method="min")

  wCalesMin2 <- calibration(data=data_employees, marginMatrix=margins, colWeights="weight"
                           , description=FALSE, method="min")

  expect_equal(wCalesMin, wCalesMin2, tolerance=1e-6)

})
