#!/bin/bash

# Exit if any command fails
set -e

# Run R script inline with correct syntax
Rscript -e '
library(infPackage)
Sys.setenv(OMP_NUM_THREADS = 4)  # Set 4 threads
system("echo $CC")          # Check the C compiler
system("echo $CXX")         # Check the C++ compiler
system("echo $CFLAGS")      # Check C flags
system("echo $CXXFLAGS")    # Check C++ flags
readFilePrint("/Users/daniel/ecological-inference-elections/instances/J100_M50_G3_I3_L50_seed18.json")
'
