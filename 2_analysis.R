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
summarize_run(res.nobias) 
# estimates: d = .00, .30, .60; group 1 Type I, 5%; Moderator power: 100%, 100%; 
# TODO: Get Type I error / power of Egger test
# PET shows not-so-bad downward bias, d = .26. Wasn't it worse before?
# check out the super poor power on the interactive model!

summarize_run(res.medPB)
# estimates: d = .13, .40, .63; group 1 Type I, 81%; Moderator power: 100%, 100%;
# note that differences between subpopulations have been halved
# PET bias is really bad, probably drug down by the bad subgroup
# interactive model makes it worse, not better.

summarize_run(res.smallfx.nobias) 
# observed d = 0, .2, .4; power = 9%, 87%, 100%; egger type 1 ??%
summarize_run(res.smallfx.medPB) # dramatically reduced power
# observed d = .14, .33, .47; power = 100%, 2.2%, 35%; egger power 74%

plotCellMeans(res.nobias, res.medPB, "No Bias", "Med Bias")
plotParamMeans(res.nobias, res.medPB, "No Bias", "Med Bias")
plotCellMeans(res.smallfx.nobias, res.smallfx.medPB, "No Bias", "Med Bias")
plotParamMeans(res.smallfx.nobias, res.smallfx.medPB, "No Bias", "Med Bias")

summarize_run(res.medPB) 
# cell means: .35, .44, .64
# Additive PET-RMA has downward bias: .00, .11, .29
# Interactive PET-RMA has downward bias, but less bad: -.10, .09, .44
# moderator power: 34%, 100%
# PET-RMA seems to make moderator power worse, not better.

summarize_run(res.hiPB)
plotCellMeans(res.nobias, res.medPB, "No Bias", "Hi Bias")
plotParamMeans(res.nobias, res.medPB, "No Bias", "Hi Bias")
# Not that bad, but these are big differences

summarize_run(res.smallfx.hiPB)
plotCellMeans(res.smallfx.nobias, res.smallfx.medPB, "No Bias", "Hi Bias")
plotParamMeans(res.smallfx.nobias, res.smallfx.medPB, "No Bias", "Hi Bias")
# power drops to 47%, 97% from 87%, 100%
# PET-RMA seems to improve power slightly (additive) but terrible in interactive

# Story not quite so damning as I'd thought if we conduct the simulations this way
# QRPs may change the story but I'm not sure if I've got the cycles
# Check http://shinyapps.org/apps/metaExplorer/ to see if things look like real-life metas
# I might say it does. Although maybe med + high QRP is more like it?
# RE with med bias, no QRP turns delta = 0, 0.2, 0.5 to d = .15, .32, .55
# RE with med bias, high QRP turns delta = 0, 0.2, 0.5 to d = .28, .38, .57
# RE with hi bias, no QRP turns delta = 0, 0.2, 0.5 to d = .3, .38, .56
# RE with hi bias hi QRP turns delta = 0, 0.2, 0.5 to d = .39, .41, .57
# It doesn't seem to matter a whole lot for the naive estimator. but QRP will mess up PET.


# Next steps:
# increase nSim
# inspect intermediate steps for quality, e.g. unit tests
# explore different values of vector d
# check distribution of N -- is it appropriate?
hist(meta1$N)
# implement PET-RMA model


