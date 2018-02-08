# Load packages and functions
library(MASS)
library(tidyverse)
library(truncnorm)
library(truncdist)
library(pwr)
library(compiler)
library(metafor)

# load study-simulating functions
source("sim-studies.R")

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
set.seed(42069)
testset <- modMA(k = 20, d = c(0, .3, .6), 
                 sel = 0, propB = 0, QRP = 0)
# check contrasts
testset$id
# make funnels for each subgroup
funnel(rma(yi = d, sei = se, data = testset, subset = id == 1))
funnel(rma(yi = d, sei = se, data = testset, subset = id == 2))
funnel(rma(yi = d, sei = se, data = testset, subset = id == 3))
# simulate biased data
testset.70 <- modMA(k = 20, d = c(0, .3, .6), 
                    sel = 1, propB = .70, QRP = 0)
# check contrasts
testset.70$id
# make funnels for each subgroup
funnel(rma(yi = d, sei = se, data = testset.70, subset = id == 1))
funnel(rma(yi = d, sei = se, data = testset.70, subset = id == 2))
funnel(rma(yi = d, sei = se, data = testset.70, subset = id == 3))

# Run simulations ----
# TODO: It would also be smart to save the output of summarize_run to an object
#       but maybe we work with the raw output objects for now

# set 1 : no bias
output.nobias <- runStudy(nSim = 100, k = 20, d = c(0, .3, .6), 
                          sel = 0, propB = 0, QRP = 0)
output.nobias

# set 2: publication bias, 70+% of each subgroup stat. sig
output.70p_pubbias <- runStudy(nSim = 100, k = 20, d = c(0, .3, .6), 
                               sel = 1, propB = .70, QRP = 0)
output.70p_pubbias

# set 3: subtler differences
output.smallfx.nobias <- runStudy(nSim = 100, k = 20, d = c(0, .2, .4), 
                                  sel = 0, propB = 0, QRP = 0)

# set 4: subtler differences + 70% of results are stat. sig
output.smallfx.70p_pubbias <- runStudy(nSim = 100, k = 20, d = c(0, .2, .4), 
                                       sel = 1, propB = 0.7, QRP = 0)

# summarize results of set 1 and set 2
summarize_run(output.nobias) 
# estimates: d = .00, .30, .60; power, 4.2%, 92%, 100%; egger type 1, 11%
# PET shows strong downward bias, d = .14

summarize_run(output.70p_pubbias)
# estimates: d = .45, .58, .74; power, 100%, 15%, 91%; egger power, 76%
# note that differences between subpopulations have been halved

summarize_run(output.smallfx.nobias) 
# observed d = 0, .2, .4; power = 5.6%, 60%, 99%; egger type 1 6.4%
summarize_run(output.smallfx.70p_pubbias) # dramatically reduced power
# observed d = .45, .54, .63; power = 100%, 2.2%, 35%; egger power 74%





# This area under construction

# Plot
output.nobias$bias <- "No Bias"
output.70p_pubbias$bias <- "70% sig"
# Parameterized as cell means (contrast coding)
bind_rows(output.nobias, output.70p_pubbias) %>% 
  mutate(d1.obs = mod.b.obs.1 + mod.b.obs.2,
         d2.obs = mod.b.obs.1 + mod.b.obs.3,
         d3.obs = mod.b.obs.1 - mod.b.obs.2 - mod.b.obs.3) %>% 
  gather(key, value, d1.obs:d3.obs) %>% 
  ggplot(aes(x = value)) +
  geom_histogram() +
  facet_grid(bias ~ key) +
  scale_x_continuous(limits = c(-.2, 1))

# Parameterized as moderator's effect size
bind_rows(output.nobias, output.70p_pubbias) %>% 
  gather(key, value, mod.b.obs.2:mod.b.obs.3) %>% 
  ggplot(aes(x = value)) +
  geom_histogram() +
  facet_grid(bias ~ key) +
  scale_x_continuous(limits = c(-1, .5))  +
  ggtitle("Moderation between d = 0, 0.3, 0.6")

# Next steps:
# increase nSim
# inspect intermediate steps for quality, e.g. unit tests
# explore different values of vector d
# check distribution of N -- is it appropriate?
hist(meta1$N)
# implement PET-RMA model

# Plot
output.smallfx.nobias$bias <- "No Bias"
output.smallfx.70p_pubbias$bias <- "70% sig"
# Parameterized as cell means
bind_rows(output.smallfx.nobias, output.smallfx.70p_pubbias) %>% 
  mutate(d1.obs = mod.b.obs.1 + mod.b.obs.2,
         d2.obs = mod.b.obs.1 + mod.b.obs.3,
         d3.obs = mod.b.obs.1 - mod.b.obs.2 - mod.b.obs.3) %>% 
  gather(key, value, d1.obs:d3.obs) %>% 
  ggplot(aes(x = value)) +
  geom_histogram() +
  facet_grid(bias ~ key) +
  scale_x_continuous(limits = c(-.2, 1))

# Parameterized as moderator's effect size
bind_rows(output.smallfx.nobias, output.smallfx.70p_pubbias) %>% 
  gather(key, value, mod.b.obs.2:mod.b.obs.3) %>% 
  ggplot(aes(x = value)) +
  geom_histogram() +
  facet_grid(bias ~ key) +
  scale_x_continuous(limits = c(-.2, 1)) +
  ggtitle("Moderation between d = 0, 0.2, 0.4")




# playing with PET-RMA
petrma.easy <- modMA(100, d = c(0, 1, 2),
                QRP = 0, sel = 0, propB = 0)
petrma.add.easy <- rma(yi = d, sei = se, 
                       mods = ~ id + se, data = petrma.easy)
summary(petrma.add.easy)
petrma.mul.easy <- rma(yi = d, sei = se, 
                       mods = ~ id * se, data = petrma.easy)
summary(petrma.mul.easy)

petrma.bias <- modMA(100, d = c(0, 1, 2),
                     QRP = 0, sel = 1, propB = .7)
petrma.add.bias <- rma(yi = d, sei = se, 
                       mods = ~ id + se, data = petrma.bias)
summary(petrma.add.bias)
petrma.mul.bias <- rma(yi = d, sei = se, 
                       mods = ~ id * se, data = petrma.bias)
summary(petrma.mul.bias)
# 
# # fitting model to k = 3000 seems very computationally taxing somehow
# petrma.easy.1k <- modMA(1000, d = c(0, 1, 2),
#                      QRP = 0, sel = 0, propB = 0)
# petrma.add.easy.1k <- rma(yi = d, sei = se, 
#                        mods = ~ id + se, data = petrma.easy.1k)
# summary(petrma.add.easy.1k) 
# # estimates d = -.25, .97, 1.9
# # significant egger term, b = .92, p < .0001. silly.
# petrma.mul.easy.1k <- rma(yi = d, sei = se, 
#                        mods = ~ id * se, data = petrma.easy.1k)
# summary(petrma.mul.easy.1k)
# # estimates d = 0, .7, 1.5
# # sig egger term for the larger effects, maybe because of d rather than g
# # interpretation of terms seems as one would expect. It's just that
# # it's very sensitive and probably rather biased and clumsy
# 
# petrma.bias.1k <- modMA(1000, d = c(0, 1, 2),
#                      QRP = 0, sel = 1, propB = .7)
# petrma.add.bias.1k <- rma(yi = d, sei = se, 
#                        mods = ~ id + se, data = petrma.bias.1k)
# summary(petrma.add.bias.1k)
# petrma.mul.bias.1k <- rma(yi = d, sei = se, 
#                        mods = ~ id * se, data = petrma.bias.1k)
# summary(petrma.mul.bias.1k)