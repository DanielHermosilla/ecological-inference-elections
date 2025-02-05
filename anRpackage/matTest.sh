#!/bin/bash

# Exit if any command fails
set -e

# Run R script inline with correct syntax
R --vanilla -e '
library(infPackage)

for (s in 1:20)
	{
		file <- paste0("/Users/daniel/ecological-inference-elections/instances/J100_M50_G3_I3_L50_seed", s , ".json")
		readFromFile(file)
		results <- EMAlgorithmAll("Multinomial")
		print(results)
	}
'
