#!/bin/bash

# Exit if any command fails
set -e

# Run R script inline with correct syntax
R --vanilla -e '
library(infPackage)

for (i in 1:21) {
	a <- EMModel$new(jsonPath = paste0("/Users/daniel/ecological-inference-elections/anRpackage/instances/J100_M50_G2_I3_L50_seed",i,".json"))
a$compute(main_method = "MVN PDF")
a$compute(main_method = "Multinomial")
a$compute(main_method = "MVN CDF")
a$compute(main_method = "Hit and Run", samples = 1000, step_size = 3000)
a$compute(main_method = "Exact")
cat(paste0("Ready with the seed: ",i))
print(class(a))
summary(a)
Sys.sleep(3)
}
'
