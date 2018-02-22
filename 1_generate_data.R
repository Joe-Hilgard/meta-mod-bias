# Generate simulated datasets

# Load packages and functions
library(MASS)
library(tidyverse)
library(truncnorm)
library(truncdist)
library(pwr)
library(compiler)
library(metafor)

# load study-simulating functions
source("sim-studies/sim-studies.R", chdir=TRUE)
source("hilgard_functions.R")

# Tough choice between dataMA and simMA.
# dataMA did the old approach of letting you force some percentage of studies to be stat sig.
# simMA does the new approach of generating a bunch of studies, then publishing only some percentage of the nulls.
# For simMA results to be noticeably biased or look anything like real-world data, you have to implement QRPs.
# or really strong selection bias. And QRPs require a lot of cycles.

# Maybe the most appropriate simulation settings are to use simMA with both QRPs and heavy pub bias.
# That's probably what's going on in the real world.

# Steps: 
# 1. simulate studies
# 2. filter for significance
# 3. run moderator model w/ and w/o filter

# To simulate studies, use dataMA()
# I'm a little surprised dataMA doesn't generate Hedges' g

# make example data set for testing & debugging
# set.seed(42069)
# meta1 <- modMA(k = 20, d = c(0, .3, .6), 
#                sel = 1, propB = .70, QRP = 0)
# m1 <- rma(yi = d, sei = se, data = meta1)
# funnel(m1, pch = as.numeric(meta1$id))
# m1.mod <- rma(yi = d, sei = se, mods = ~id, data = meta1)
# summary(m1.mod)
# funnel(m1.mod)
# m1.egger <- rma(yi = d, sei = se, mods = ~se, data = meta1)
# summary(m1.egger)

# For starters, we can ask:
# By how much do we miss mean effect size of mean(d %*% k)
# What is the power to detect d ~ id?
# By how much does that power drop when we filter for sig?

# testing that passing of arguments is working right
set.seed(999)
testset <- modMA(k = 30, d = c(0, .3, .6))
# check contrasts
testset$id
# make funnels for each subgroup
myFunnel(rma(yi = d, sei = se, data = testset, subset = id == 1))
myFunnel(rma(yi = d, sei = se, data = testset, subset = id == 2))
myFunnel(rma(yi = d, sei = se, data = testset, subset = id == 3))
# simulate biased data
testset.med <- modMA(k = 30, d = c(0, .3, .6), 
                     censor = "med")
# check contrasts
testset.med$id
# make funnels for each subgroup
myFunnel(rma(yi = d, sei = se, data = testset.med, subset = id == 1))
myFunnel(rma(yi = d, sei = se, data = testset.med, subset = id == 2))
myFunnel(rma(yi = d, sei = se, data = testset.med, subset = id == 3))

# Run simulations ----
# TODO: It would also be smart to save the output of summarize_run to an object
#       but maybe we work with the raw output objects for now

# set 1 : no bias
output.nobias <- runStudy(nSim = 100, k = 20, d = c(0, .3, .6))
output.nobias

# set 2: publication bias, medium censoring
output.med_pubbias <- runStudy(nSim = 100, k = 20, d = c(0, .3, .6), 
                               censor = "med")
output.med_pubbias # oh boy, is convergence gonna be a thing now?

# set 3: subtler differences
output.smallfx.nobias <- runStudy(nSim = 100, k = 20, d = c(0, .2, .4))

# set 4: subtler differences + 70% of results are stat. sig
output.smallfx.med_pubbias <- runStudy(nSim = 100, k = 20, d = c(0, .2, .4), 
                                       censor = "med")

# set 5: Heavier pub bias and some QRPs.
# oh my god this takes forever to run. Simulating 20*3*100 p-hacked studies...!
# should probably generate all the data objects in one script and save them as .RData
output.bias <- runStudy(nSim = 100, k = 20, d = c(0, .3, .6), 
                        censor = "high", qrpEnv = "medium")

save.image("sim.RData")