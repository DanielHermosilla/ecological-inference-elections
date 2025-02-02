#!/bin/bash

# Exit if any command fails
set -e

# Run R script inline with correct syntax
Rscript -e '
library(infPackage)
readFilePrint("/Users/daniel/ecological-inference-elections/instances/J100_M50_G3_I3_L50_seed18.json")
'
