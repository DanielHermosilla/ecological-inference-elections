#!/bin/bash

# Exit if any command fails
set -e

# Run R script inline with correct syntax
R --vanilla -e '
library(infPackage)

a <- EMModel$new(jsonPath = "/Users/daniel/ecological-inference-elections/anRpackage/instances/J100_M50_G2_I3_L50_seed18.json")
a$compute(main_method = "MVN PDF")
Sys.sleep(3)
a$compute(main_method = "Multinomial")
Sys.sleep(3)
a$compute(main_method = "MVN CDF")
Sys.sleep(3)
a$compute(main_method = "Hit and Run", samples = 1000, step_size = 3000)
Sys.sleep(3)
a$compute(main_method = "Exact")

'
