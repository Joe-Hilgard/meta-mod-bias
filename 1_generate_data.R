# Generate simulated datasets

# Load packages and functions
library(MASS)
library(tidyverse)
library(truncnorm)
library(truncdist)
library(pwr)
library(compiler)
library(metafor)
library(weightr)

# load study-simulating functions
source("sim-studies/sim-studies.R", chdir=TRUE)
source("hilgard_functions.R")
source("tidy_functions.R")

# Tough choice between dataMA and simMA.
# dataMA did the old approach of letting you force some percentage of studies to be stat sig.
# simMA does the new approach of generating a bunch of studies, then publishing only some percentage of the nulls.
# For simMA results to be noticeably biased or look anything like real-world data, you have to implement QRPs.
# or really strong selection bias. And QRPs require a lot of cycles.

# Maybe the most appropriate simulation settings are to use simMA with both QRPs and heavy pub bias.
# That's probably what's going on in the real world.
# Do I have anything like enough processing power to do that?

# could also take data from meta-showdown and stitch together biased studies from various delta
# But that would require the raw data, which is too large to fit into RAM.

# I wonder if/how I should implement meta-moderation by a continuous covariate.

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
# check summary stats of sample size
summary(testset$N) 
# seems a bit high for social psych. also max N is over 200? Are arguments passing correctly?
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
res.nobias <- runStudy(nSim = 1000, k = 20, d = c(0, .3, .6))
res.nobias

# set 2: publication bias, medium censoring
res.medPB <- runStudy(nSim = 1000, k = 20, d = c(0, .3, .6), 
                               censor = "med")
res.medPB # oh boy, is convergence gonna be a thing now?

# set 3: subtler differences
res.smallfx.nobias <- runStudy(nSim = 100, k = 20, d = c(0, .2, .4))

# set 4: subtler differences + 70% of results are stat. sig
res.smallfx.medPB <- runStudy(nSim = 100, k = 20, d = c(0, .2, .4), 
                                       censor = "med")

# set 5: Heavy pub bias.
# runs nice and quickly. It's simulating the QRPs that gets you.
# having surprising convergence issues with plain PET
res.hiPB <- runStudy(nSim = 1000, k = 20, d = c(0, .3, .6),
                     censor = "high")

# set 6: subtler differences + heavy pub bias
res.smallfx.hiPB <- runStudy(nSim = 1000, k = 20, d = c(0, .2, .4),
                             censor = "high")

# set 7: Med pub bias and hi QRPs
# oh my god this takes forever to run. Simulating 20*3*100 p-hacked studies...!
# should probably generate all the data objects in one script and save them as .RData
res.medPB.hiQRP <- runStudy(nSim = 1000, k = 20, d = c(0, .3, .6), 
                         censor = "medium", qrpEnv = "high")

# set 8: subtler effects + med pub bias + hi QRPs
res.smallfx.medPB.hiQRP <- runStudy(nSim = 1000, k = 20, d = c(0, .2, .4), 
                                    censor = "medium", qrpEnv = "high")

# One-time kludge for fixing column names
# res.medPB.hiQRP <- res.medPB.hiQRP %>% 
#   rename(mod.obs.b1 = mod.b.obs.1,
#          mod.obs.b2 = mod.b.obs.2,
#          mod.obs.b3 = mod.b.obs.3,
#          mod.p1 = mod.p.1,
#          mod.p2 = mod.p.2,
#          mod.p3 = mod.p.3)
# 
# res.nobias <- res.nobias %>% 
#   rename(mod.obs.b1 = mod.b.obs.1,
#          mod.obs.b2 = mod.b.obs.2,
#          mod.obs.b3 = mod.b.obs.3,
#          mod.p1 = mod.p.1,
#          mod.p2 = mod.p.2,
#          mod.p3 = mod.p.3)
# 
# res.smallfx.nobias <- res.smallfx.nobias %>% 
#   rename(mod.obs.b1 = mod.b.obs.1,
#          mod.obs.b2 = mod.b.obs.2,
#          mod.obs.b3 = mod.b.obs.3,
#          mod.p1 = mod.p.1,
#          mod.p2 = mod.p.2,
#          mod.p3 = mod.p.3)
# 
# res.smallfx.medPB.hiQRP <- res.smallfx.medPB.hiQRP %>% 
#   rename(mod.obs.b1 = mod.b.obs.1,
#          mod.obs.b2 = mod.b.obs.2,
#          mod.obs.b3 = mod.b.obs.3,
#          mod.p1 = mod.p.1,
#          mod.p2 = mod.p.2,
#          mod.p3 = mod.p.3)

save.image("sim.RData")

# TODO: Implement selection modeling
# Without moderator
with(testset,
     weightfunct(effect = d, v = v))
# With moderator
with(testset,
     weightfunct(effect = d, v = v, mods = ~id))
