# copyright (C) 2014-2016 A.Rebecq
library(testthat)

context("Test calibration functions on big dataset")

test_that("Calibration functions check out with Calmar", {
  
  data("population_test")
  
  sample <- dataPop[dataPop$weight > 0,]

  wCalLin <- calibration(data=sample, marginMatrix=table_margins_1, colWeights="weight"
                         , method="linear", description=FALSE)
  
  ## Just checking that calibration works with linear method
  ## when all variables are categorical
  wCalLin_categ <- calibration(data=sample, marginMatrix=table_margins_1[10:11,], colWeights="weight"
                               , method="linear", description=FALSE)
  
  expect_length(wCalLin_categ, 1000)
  
  wCalRaking <- calibration(data=sample, marginMatrix=table_margins_1, colWeights="weight"
                         , method="raking", description=FALSE, pct=FALSE)
  
  wCalLogit <- calibration(data=sample, marginMatrix=table_margins_1, colWeights="weight"
                            , method="logit", bounds=c(0.2,1.3), description=FALSE)
  
  
  expect_equal(wCalLin, poptest_calmar$weight_cal_lin, tolerance=1e-6)
  expect_equal(wCalRaking, poptest_calmar$weight_cal_raking, tolerance=1e-6)
  expect_equal(wCalLogit, poptest_calmar$weight_cal_logit, tolerance=1e-6)
  
  popTotal <- 50000
  
  wCalLin2 <- calibration(data=sample, marginMatrix=table_margins_2, colWeights="weight"
                         , method="linear", description=FALSE, popTotal=popTotal, pct=TRUE)
  
  wCalRaking2 <- calibration(data=sample, marginMatrix=table_margins_2, colWeights="weight"
                            , method="raking", description=FALSE, popTotal=popTotal, pct=TRUE)
  
  wCalLogit2 <- calibration(data=sample, marginMatrix=table_margins_2, colWeights="weight"
                           , method="logit", bounds=c(0.2,1.3), description=FALSE, popTotal=popTotal, pct=TRUE)
  
  
  expect_equal(wCalLin2, poptest_calmar$weight_cal_lin_2, tolerance=1e-6)
  expect_equal(wCalRaking2, poptest_calmar$weight_cal_raking_2, tolerance=1e-6)
  expect_equal(wCalLogit2, poptest_calmar$weight_cal_logit_2, tolerance=1e-6)
  
  ## Test that calibrated estimators are equal to population totals
  totalsCalVar <- sapply(table_margins_1[,1], function(x) { return(sum(dataPop[,x])) })
  expect_equal( sapply(table_margins_1[,1], function(x) { return(HTtotal(sample[,x], wCalLin)) }),
                totalsCalVar,
                tolerance=1e-6)
  expect_equal( sapply(table_margins_1[,1], function(x) { return(HTtotal(sample[,x], wCalLogit)) }),
                totalsCalVar,
                tolerance=1e-6)
  expect_equal( sapply(table_margins_1[,1], function(x) { return(HTtotal(sample[,x], wCalRaking2)) }),
                totalsCalVar,
                tolerance=1e-6)
  
  ########### Tests with non-response and scale factor
  ## TODO : check parameters for logit calibration
  sampleNR <- sample[sample$responding==1,]
  
  wCalLinNR <- calibration(data=sampleNR, marginMatrix=table_margins_1, colWeights="weight"
                         , method="linear", description=FALSE, scale=TRUE)
  wCalRakingNR <- calibration(data=sampleNR, marginMatrix=table_margins_1, colWeights="weight"
                           , method="raking", description=FALSE, scale=TRUE)
  wCalLogitNR <- calibration(data=sampleNR, marginMatrix=table_margins_1, colWeights="weight"
                           , method="logit", bounds=c(0.1,2.0), description=FALSE, scale=TRUE, popTot=50000)
  
  expect_equal(wCalLinNR, poptest_calmar_nr$weight_cal_lin, tolerance=1e-6)
  expect_equal(wCalRakingNR, poptest_calmar_nr$weight_cal_raking, tolerance=1e-6)
  expect_equal(wCalLogitNR, poptest_calmar_nr$weight_cal_logit, tolerance=1e-6)
  
#   testDistrib <- poptest_calmar_nr$weight_cal_logit/sampleNR$weight * sum(sampleNR$weight) / 50000 
#   print(summary(testDistrib))
  
  wCalLinNR2 <- calibration(data=sampleNR, marginMatrix=table_margins_2, popTot = 50000, pct=T, colWeights="weight"
                           , method="linear", description=FALSE, scale=TRUE)
  wCalRakingNR2 <- calibration(data=sampleNR, marginMatrix=table_margins_2, popTot = 50000, pct=T, colWeights="weight"
                              , method="raking", description=FALSE, scale=TRUE)
  wCalLogitNR2 <- calibration(data=sampleNR, marginMatrix=table_margins_2, popTot = 50000, pct=T, colWeights="weight"
                             , method="logit", bounds=c(0.1,2.0), description=FALSE, scale=TRUE)
  
  expect_equal(wCalLinNR2, poptest_calmar_nr$weight_cal_lin_2, tolerance=1e-6)
  expect_equal(wCalRakingNR2, poptest_calmar_nr$weight_cal_raking_2, tolerance=1e-6)
  expect_equal(wCalLogitNR2, poptest_calmar_nr$weight_cal_logit_2, tolerance=1e-6)
  
  ## Checks for qk vectors
  wCalLin_q1 <- calibration(data=sample, marginMatrix=table_margins_1, colWeights="weight"
                         , q=rep(1,nrow(sample)), method="linear", description=FALSE)
  
  wCalLinNR_q1 <- calibration(data=sampleNR, marginMatrix=table_margins_1, colWeights="weight"
                           , q=rep(1,nrow(sampleNR)), method="linear", description=FALSE, scale=TRUE)
  
  expect_equal(wCalLin_q1, poptest_calmar$weight_cal_lin, tolerance=1e-6)
  expect_equal(wCalLinNR_q1, poptest_calmar_nr$weight_cal_lin, tolerance=1e-6)
  
  wCalLin_qTest <- calibration(data=sample, marginMatrix=table_margins_1, colWeights="weight"
                            , q=sample$qTest, method="linear", description=FALSE)
  
  wCalRaking_qTest <- calibration(data=sampleNR, marginMatrix=table_margins_2, colWeights="weight",
                               q=sampleNR$qTest, method="raking", description=FALSE, 
                               scale = TRUE, pct=TRUE, popTotal = 50000)
  
  expect_equal(wCalLin_qTest, poptest_calmar$weight_cal_lin_qtest, tolerance=1e-2)
  expect_equal(wCalRaking_qTest, poptest_calmar_nr$weight_cal_raking_2_qtest, tolerance=1e-2)
  
  ## Check that warning is correctly thrown when
  ## user selects an incorrect method
  expect_warning(
    calibration(data=sample, marginMatrix=table_margins_1, colWeights="weight"
                                     , method="truncated", description=FALSE, popTotal = 50000),
    "not implemented", ignore.case=T)
  
  expect_warning(
    calibration(data=sample, marginMatrix=table_margins_1, colWeights="weight"
                , method="randomstuff", description=FALSE, popTotal = 50000),
    "not implemented", ignore.case=T)
  
  expect_warning(
    calibration(data=sample, marginMatrix=table_margins_1, colWeights="weight"
                , method=NULL, description=FALSE, popTotal = 50000),
    "not specified", ignore.case=T)
  
  ## Check that errors are correctly thrown when impossible
  ## calibration margins are entered
  expect_error(
    calibration(data=sample, marginMatrix=table_margins_1, colWeights="weight"
                , method="linear", description=FALSE, popTotal = 51000, maxIter=100),
    "no convergence", ignore.case=T)
  
  expect_error(
    calibration(data=sample, marginMatrix=table_margins_1, colWeights="weight"
                , method="raking", description=FALSE, popTotal = 51000, maxIter=100),
    "no convergence", ignore.case=T)
  
  expect_error(
    calibration(data=sample, marginMatrix=table_margins_1, colWeights="weight"
                , method="logit", bounds=c(0.5,1.5), description=FALSE, 
                popTotal = 51000, maxIter=100),
    "no convergence", ignore.case=T)
  
  ## Check that bounds are correctly required for logit
  expect_error(
    calibration(data=sample, marginMatrix=table_margins_1, colWeights="weight"
              , method="logit", description=FALSE),
    "must enter LO and UP bounds", ignore.case=T)
  
  ## Check errors when q is not valid
  expect_error(
    calibration(data=sample, marginMatrix=table_margins_1, colWeights="weight"
                , q=rep(1,13), method="linear", description=FALSE),
    "Vector q must have same length as data", ignore.case=T)
  
  expect_error(
    calibration(data=sample, marginMatrix=table_margins_1, colWeights="weight"
                , q=sample$qTest, method="min", description=FALSE),
    "not supported", ignore.case=T)
  
  expect_error(
    calibration(data=sample, marginMatrix=table_margins_1, colWeights="weight"
                , q=sample$qTest, costs = rep(1,nrow(sample)), description=FALSE),
    "not supported", ignore.case=T)
  
  ## Test when margins of categorical variables are
  ## entered in percentages whose sum is 100 instead of 1
  table_margins_3 <- table_margins_2
  table_margins_3[10,3:7] <- c(20,20,20,20,20)
  table_margins_3[11,3:5] <- c(10,60,30)
  
  wCalLin3 <- calibration(data=sample, marginMatrix=table_margins_3, colWeights="weight"
                          , method="linear", description=FALSE, popTotal=popTotal, pct=TRUE)
  
  wCalRaking3 <- calibration(data=sample, marginMatrix=table_margins_3, colWeights="weight"
                             , method="raking", description=FALSE, popTotal=popTotal, pct=TRUE)
  
  wCalLogit3 <- calibration(data=sample, marginMatrix=table_margins_3, colWeights="weight"
                            , method="logit", bounds=c(0.2,1.3), description=FALSE, popTotal=popTotal, pct=TRUE)
  
  expect_equal(wCalLin2, wCalLin3, tolerance=1e-6)
  expect_equal(wCalRaking2, wCalRaking3, tolerance=1e-6)
  expect_equal(wCalLogit2, wCalLogit3, tolerance=1e-6)
  
})

test_that("Test margin stats", {
  
  data("population_test")
  
  sample <- dataPop[dataPop$weight > 0,]
  
  testStats1 <- marginStats(sample, table_margins_1, colWeights = "weight")
  testStats2 <- marginStats(sample, table_margins_2, colWeights = "weight", pct=T, popTotal = 50000)
  
  expect_equal(testStats1[13,1], 18.94)
  expect_equal(testStats1[14,2], 20.02)
  expect_equal(testStats1[15,3], 1.44)
  
  expect_equal(testStats2[17,3], 0.27)
  expect_equal(testStats2[18,1], 30.21)
  expect_equal(testStats2[16,2], 10)
  
  ## TODO : test marginStats with calibration weights
  ## (on penalized calibration for instance)
  sample$wCal <- calibration(data=sample, marginMatrix=table_margins_1, colWeights="weight"
                      , method="linear", description=FALSE, popTotal = 50000)
  
  testCosts <- rep(Inf, length(table_margins_1[,1]))

  expect_warning(
    sample$wCal_penal <- calibration(data=sample, marginMatrix=table_margins_1, colWeights="weight"
                             , method="linear", description=FALSE, costs=testCosts, popTotal = 50000),
    "all costs are infinite", ignore.case=T)
  
  testStats3 <- marginStats(sample, table_margins_1, colWeights = "weight", colCalibratedWeights = "wCal", popTotal = 50000)
  expect_equal(testStats3[,4], rep(0, length(testStats3[,1])))
  
  testStats4 <- marginStats(sample, table_margins_1, colWeights = "weight", colCalibratedWeights = "wCal_penal", popTotal = 50000)
  expect_equal(testStats4[,4], rep(0, length(testStats4[,1])))
  
  #### TODO : improve testing for penalized calib
#   testCosts2 <- rep(Inf, length(table_margins_1[,1]))
#   testCosts2[1] <- 1
#   testCosts2[2] <- 1
#   testCosts2[3] <- 100
#   testCosts2[4] <- 1
#   testCosts2[5] <- 1
#   testCosts2[2] <- 100
#   testCosts2[3] <- 1
#   
#   print(testCosts2)
  
  # TODO : set gap
#   sample$wCal_penal2 <- calibration(data=sample, marginMatrix=table_margins_1[c(1:3),], colWeights="weight"
#                                    , method="linear", description=TRUE, costs=testCosts2[c(1:3)], popTotal=50000, gap=0.1, uCostPenalized = 1e-3)
# 
#   sample$wCal_penal2 <- calibration(data=sample, marginMatrix=table_margins_1, colWeights="weight"
#                                     , method="linear", description=TRUE, costs=testCosts2, popTotal=50000, gap=0.2, uCostPenalized = 1e-4)
#   
  
  # testStats5 <- marginStats(sample, table_margins_1[c(1:3),], colWeights = "weight", colCalibratedWeights = "wCal_penal2", popTotal = 50000)
  # testStats5 <- marginStats(sample, table_margins_1, colWeights = "weight", colCalibratedWeights = "wCal_penal2", popTotal = 50000)

  # print(testStats5)
  
})