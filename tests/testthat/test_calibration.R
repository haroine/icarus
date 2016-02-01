# copyright (C) 2014-2016 A.Rebecq
library(testthat)

context("Test calibration functions on big dataset")

test_that("Calibration functions check out with Calmar", {
  
  data("population_test")
  
  sample <- dataPop[dataPop$weight > 0,]

  wCalLin <- calibration(data=sample, marginMatrix=table_margins_1, colWeights="weight"
                         , method="linear", description=FALSE)
  
  wCalRaking <- calibration(data=sample, marginMatrix=table_margins_1, colWeights="weight"
                         , method="raking", description=FALSE)
  
  wCalLogit <- calibration(data=sample, marginMatrix=table_margins_1, colWeights="weight"
                            , method="logit", bounds=c(0.2,1.3), description=FALSE)
  
  
  expect_equal(wCalLin, poptest_calmar$weight_cal_lin, tolerance=1e-6)
  expect_equal(wCalRaking, poptest_calmar$weight_cal_raking, tolerance=1e-6)
  expect_equal(wCalLogit, poptest_calmar$weight_cal_logit, tolerance=1e-6)
  
  popTotal <- 50000
  
  wCalLin2 <- calibration(data=sample, marginMatrix=table_margins_2, colWeights="weight"
                         , method="linear", description=FALSE, popTotal=popTotal)
  
  wCalRaking2 <- calibration(data=sample, marginMatrix=table_margins_2, colWeights="weight"
                            , method="raking", description=FALSE, popTotal=popTotal)
  
  wCalLogit2 <- calibration(data=sample, marginMatrix=table_margins_2, colWeights="weight"
                           , method="logit", bounds=c(0.2,1.3), description=FALSE, popTotal=popTotal)
  
  
  expect_equal(wCalLin2, poptest_calmar$weight_cal_lin_2, tolerance=1e-6)
  expect_equal(wCalRaking2, poptest_calmar$weight_cal_raking_2, tolerance=1e-6)
  expect_equal(wCalLogit2, poptest_calmar$weight_cal_logit_2, tolerance=1e-6)
  
})
