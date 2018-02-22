# Analyse results

# Load packages and functions
library(MASS)
library(tidyverse)
library(truncnorm)
library(truncdist)
library(pwr)
library(compiler)
library(metafor)
source("sim-studies/sim-studies.R", chdir=TRUE)
source("hilgard_functions.R")

# Load simulated data
load("sim.Rdata")

# summarize results of set 1 and set 2
summarize_run(output.nobias) 
# estimates: d = .00, .30, .60; group 1 Type I, 5%; Moderator power: 100%, 100%; 
# TODO: Get Type I error / power of Egger test
# PET shows not-so-bad downward bias, d = .26. Wasn't it worse before?
# check out the super poor power on the interactive model!

summarize_run(output.med_pubbias)
# estimates: d = .13, .40, .63; group 1 Type I, 81%; Moderator power: 100%, 100%;
# note that differences between subpopulations have been halved
# PET bias is really bad, probably drug down by the bad subgroup
# interactive model makes it worse, not better.

summarize_run(output.smallfx.nobias) 
# observed d = 0, .2, .4; power = 9%, 87%, 100%; egger type 1 ??%
summarize_run(output.smallfx.med_pubbias) # dramatically reduced power
# observed d = .14, .33, .47; power = 100%, 2.2%, 35%; egger power 74%

plotCellMeans(output.nobias, output.med_pubbias, "No Bias", "Med Bias")
plotParamMeans(output.nobias, output.med_pubbias, "No Bias", "Med Bias")
plotCellMeans(output.smallfx.nobias, output.smallfx.med_pubbias, "No Bias", "Med Bias")
plotParamMeans(output.smallfx.nobias, output.smallfx.med_pubbias, "No Bias", "Med Bias")

summarize_run(output.bias) 
# cell means: .35, .44, .64
# Additive PET-RMA has downward bias: .00, .11, .29
# Interactive PET-RMA has downward bias, but less bad: -.10, .09, .44
# moderator power: 34%, 100%
# PET-RMA seems to make moderator power worse, not better.

# Next steps:
# increase nSim
# inspect intermediate steps for quality, e.g. unit tests
# explore different values of vector d
# check distribution of N -- is it appropriate?
hist(meta1$N)
# implement PET-RMA model


