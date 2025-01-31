#!/bin/bash

# Exit if any command fails
set -e

# Run R script inline with correct syntax
Rscript -e '
library(infPackage)
matrix_1 <- matrix(c(
  1.0, 3.0, 5.0,
  2.0, 4.0, 6.0,
  7.0, 8.0, 9.0
), nrow = 3, byrow = TRUE)

matrix_2 <- matrix(as.double(1:24), nrow = 3, ncol = 8)

RsetParameters(matrix_1, matrix_2)
'

