# Standard testthat entry point for R CMD check.
# Loads the installed package and runs the test suite under tests/testthat/.
library(testthat)
library(HLAtools)

test_check("HLAtools")
