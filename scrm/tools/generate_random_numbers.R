#!/usr/bin/Rscript --vanilla
#
# generate_random_numbers
# Generates a file with random numbers drawn from a 
# uniform distribution over [0, 1]
# 
# Author:   Paul R. Staab 
# Email:    staab (at) bio.lmu.de
# Date:     2012-08-30
# Licence:  GPLv3 or later
#

args <- commandArgs(TRUE)

# Enable just-in-time-compiler if availible
if ("compiler" %in% rownames(installed.packages())){
  library("compiler")
  invisible(compiler::enableJIT(3))
}

n <- 2^20
random <- runif(n)
random <- matrix(random, length(random), 1)
write.table(random, file="random.numbers", row.names=F, col.names=F)
