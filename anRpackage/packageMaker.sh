#!/bin/bash

# Exit if any command fails
set -e

# Remove old compiled files
rm -f *.tar.gz
rm -f src/RcppExports.cpp
rm -f R/RcppExports.R

# Load Rcpp and run compileAttributes
Rscript -e "library(Rcpp); Rcpp::compileAttributes()"

# Build and install the package
R CMD BUILD .
R CMD INSTALL infPackage_0.1.0.tar.gz

echo "Package successfully built and installed!"
echo "The functions created are as follow:"
Rscript -e "library(infPackage); ls("package:infPackage")"

