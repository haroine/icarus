library(testthat)
#library(gaston)
# source("../../R/ginv.R", encoding = "latin1")

context("Test matrix inversion")

test_that("Matrix inversion gives acceptable result", {
 
  mata <- matrix(c(2,1,1,0), nrow = 2, ncol = 2, byrow = T)
  matb <- matrix(c(0,1,1,-2), nrow = 2, ncol = 2, byrow = T)
  
  # TODO : replace by gaston::ginv to avoid any confusion
  expect_equal(ginv(mata), matb)
  
})