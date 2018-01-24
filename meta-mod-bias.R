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

# dataMA() documentation:
#Produces a dataset for meta-analysis. Applies both QRP
#and selection at a proportion specified by propB if 
#sel and QRP are 1 not 0. 

#' @param k the number of studies in the MA
#' @param QRP 1 if QRP/p-hacks are available, 0 otherwise
#' @param sel 1 if publication bias selection exists, 0 otherwise
#' @param propB the proportion of the sample affected by bias
#' @param meanD the true effect (or the average of the true effects if heterogeneity exists)
#' @param sigma the SD around the true effect
#' @param cbdv the correlation between the multiple DVs
#' @param maxN the max possible group size that could be created *this needs to be set higher than what can actually be generated--it doesn't mean you get bigger samples
#' @param minN the min of the truncated normal for sample size
#' @param meanN the mean of the truncated normal for sample size
#' @param sdN the SD of the truncated normal for sample size
#' @param multDV 1 if multiple DVs as a hack, 0 otherwise
#' @param out 1 if optional outlier removal as a hack, 0 otherwise
#' @param mod 1 if optional moderator as a hack, 0 otherwise
#' @param colLim number of times to try collecting more data
#' @param add number to add to each group when collecting more data
#' @param verbose Should informations be printed?
#' 

# Make a function that runs dataMA three times,
# once for each population,
# then mixes them together
# TODO: Note that other parameters might be worth exploring
#       Expand argument assignment and passing
modMA <- function(k, d, QRP, sel, propB) {
  d1 <- dataMA(k = k, 
               QRP = QRP, sel = sel, propB = propB, 
               meanD = d[1], sigma = 0,
               cbdv = .5, maxN = 200, minN = 20,
               meanN = 50, sdN = 20,
               multDV = 0, out = 0, mod = 0,
               colLim = 0, add = 0, verbose = T)
  d2 <- dataMA(k = k, 
               QRP = QRP, sel = sel, propB = propB, 
               meanD = d[2], sigma = 0,
               cbdv = .5, maxN = 200, minN = 20,
               meanN = 50, sdN = 20,
               multDV = 0, out = 0, mod = 0,
               colLim = 0, add = 0, verbose = T)
  d3 <- dataMA(k = k, 
               QRP = QRP, sel = sel, propB = propB, 
               meanD = d[3], sigma = 0,
               cbdv = .5, maxN = 200, minN = 20,
               meanN = 50, sdN = 20,
               multDV = 0, out = 0, mod = 0,
               colLim = 0, add = 0, verbose = T)
  return(bind_rows(d1, d2, d3, .id = "id"))
}

# make example data set for testing & debugging
set.seed(42069)
meta1 <- modMA(k = 20, d = c(0, .3, .6), 
               sel = 1, propB = .70, QRP = 0)

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

# function to analyze results of dataMA
inspectMA <- function(dataset) {
  # basic model
  rmamod <- rma(yi = d, sei = se, data = dataset)
  # test for moderator
  moderation.test <- rma(yi = d, sei = se, 
                         mods = ~id, data = dataset)
  # test for small-study effect
  egger.test <- rma(yi = d, sei = se,
                    mods = ~se, data = dataset)
  # test for moderator after adjustment for small-study
  # is it id + se or id * se?
  joint.test.additive <- rma(yi = d, sei = se,
                    mods = ~id + se, data = dataset)
  joint.test.interactive <- rma(yi = d, sei = se,
                             mods = ~id * se, data = dataset)
  # not sure how to interpret these parameters in the interactive model
  # test for moderator in hedges & vevea weight model
  # TODO  
  # return test results
  out <- data.frame(d.obs = summary(rmamod)$b[1],
                    se.obs = summary(rmamod)$se[1],
                    d.p = summary(rmamod)$pval[1],
                    # Moderators
                    # Note that .1 is the intercept, reference group
                    mod.b.obs.1 = summary(moderation.test)$b[1],
                    mod.b.obs.2 = summary(moderation.test)$b[2],
                    mod.b.obs.3 = summary(moderation.test)$b[3],
                    mod.p.1 = summary(moderation.test)$pval[1],
                    mod.p.2 = summary(moderation.test)$pval[2],
                    mod.p.3 = summary(moderation.test)$pval[3],
                    # Egger
                    d.obs.pet = summary(egger.test)$b[1],
                    p.pet = summary(egger.test)$pval[1],
                    b.egger = summary(egger.test)$b[2],
                    p.egger = summary(egger.test)$pval[2]
                    # joint PET-RMA tests
                    )
  # to be continued...
  
  
  return(data.frame(out))
}


# testing that passing of arguments is working right
set.seed(42069)
testset <- modMA(k = 20, d = c(0, .3, .6), 
                 sel = 0, propB = 0, QRP = 0)
funnel(rma(yi = d, sei = se, data = testset, subset = id == 1))
funnel(rma(yi = d, sei = se, data = testset, subset = id == 2))
funnel(rma(yi = d, sei = se, data = testset, subset = id == 3))
testset.70 <- modMA(k = 20, d = c(0, .3, .6), 
                    sel = 1, propB = .70, QRP = 0)
funnel(rma(yi = d, sei = se, data = testset.70, subset = id == 1))
funnel(rma(yi = d, sei = se, data = testset.70, subset = id == 2))
funnel(rma(yi = d, sei = se, data = testset.70, subset = id == 3))

# Run simulations ----
nSim <- 500

# set 1 : no bias
output.nobias <- data.frame(NA); 
for (i in 1:nSim) {
  # this iteration, run dataMA() on three populations
  # k studies each, effect sizes given by d
  t <- modMA(k = 20, d = c(0, .3, .6), 
             sel = 0, propB = 0, QRP = 0) %>% 
    inspectMA()
  # add this iteration's results to the output object
  output.nobias <- bind_rows(output.nobias, t)
}
# Trim missings from first row, first column.
# Surely there's a more elegant way
output.nobias <- output.nobias[-1,]
output.nobias <- output.nobias[,-1]

# set 2: publication bias, 70+% of each subgroup stat. sig
output.70p_pubbias <- data.frame(NA)
for (i in 1:nSim) {
  # this iteration, run dataMA() on three populations
  # k studies each, effect sizes given by d
  t <- modMA(k = 20, d = c(0, .3, .6),
             sel = 1, QRP = 0, propB = .70) %>% 
    inspectMA()
  # add this iteration's results to the output object
  output.70p_pubbias <- bind_rows(output.70p_pubbias, t)
}

output.70p_pubbias
# Trim missings from first row, first column.
# Surely there's a more elegant way
output.70p_pubbias <- output.70p_pubbias[-1,]
output.70p_pubbias <- output.70p_pubbias[,-1]

# ideal case stats:
is_sig <- function(x) x < .05

summarize_run <- function(x) {
  x.est <- x %>% 
    summarize_each(funs(mean), d.obs, se.obs, mod.b.obs.1:mod.b.obs.3,
                   d.obs.pet, b.egger) %>% 
    mutate(d1.obs = mod.b.obs.1,
           d2.obs = mod.b.obs.1 + mod.b.obs.2,
           d3.obs = mod.b.obs.1 + mod.b.obs.3)

  x.pow <- x %>% 
    summarize_each(funs(mean(is_sig(.))), 
                   d.p, mod.p.1:mod.p.3, p.egger)
  
  return(bind_cols(x.est, x.pow))
}

summarize_run(output.nobias) 

# estimates: d = .00, .30, .60; power, 4.2%, 92%, 100%; egger type 1, 11%
# PET shows strong downward bias, d = .14
summarize_run(output.70p_pubbias)
# estimates: d = .45, .58, .74; power, 100%, 15%, 91%; egger power, 76%
# note that differences between subpopulations have been halved

# Plot
output.nobias$bias <- "No Bias"
output.70p_pubbias$bias <- "70% sig"
# Parameterized as cell means
bind_rows(output.nobias, output.70p_pubbias) %>% 
  mutate(d1.obs = mod.b.obs.1,
         d2.obs = mod.b.obs.1 + mod.b.obs.2,
         d3.obs = mod.b.obs.1 + mod.b.obs.3) %>% 
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
  scale_x_continuous(limits = c(-.2, 1))  +
  ggtitle("Moderation between d = 0, 0.3, 0.6")

# Next steps:
# increase nSim
# inspect intermediate steps for quality, e.g. unit tests
# explore different values of vector d
# check distribution of N -- is it appropriate?
hist(meta1$N)
# implement PET-RMA model

# what if differences are more subtle?
output.smallfx.nobias <- data.frame(NA)
for (i in 1:nSim) {
  # this iteration, run dataMA() on three populations
  # k studies each, effect sizes given by d
  t <- modMA(k = 20, d = c(0, .2, .4),
             sel = 0, QRP = 0, propB = 0) %>% 
    inspectMA()
  # add this iteration's results to the output object
  output.smallfx.nobias <- bind_rows(output.smallfx.nobias, t)
}
output.smallfx.nobias <- output.smallfx.nobias[-1,]
output.smallfx.nobias <- output.smallfx.nobias[,-1]

output.smallfx.70p_pubbias <- data.frame(NA)
for (i in 1:nSim) {
  # this iteration, run dataMA() on three populations
  # k studies each, effect sizes given by d
  t <- modMA(k = 20, d = c(0, .2, .4),
             sel = 1, QRP = 0, propB = .70) %>% 
    inspectMA()
  # add this iteration's results to the output object
  output.smallfx.70p_pubbias <- bind_rows(output.smallfx.70p_pubbias, t)
}
output.smallfx.70p_pubbias <- output.smallfx.70p_pubbias[-1,]
output.smallfx.70p_pubbias <- output.smallfx.70p_pubbias[,-1]

summarize_run(output.smallfx.nobias) 
# observed d = 0, .2, .4; power = 5.6%, 60%, 99%; egger type 1 6.4%
summarize_run(output.smallfx.70p_pubbias) # dramatically reduced power
# observed d = .45, .54, .63; power = 100%, 2.2%, 35%; egger power 74%

# Plot
output.smallfx.nobias$bias <- "No Bias"
output.smallfx.70p_pubbias$bias <- "70% sig"
# Parameterized as cell means
bind_rows(output.smallfx.nobias, output.smallfx.70p_pubbias) %>% 
  mutate(d1.obs = mod.b.obs.1,
         d2.obs = mod.b.obs.1 + mod.b.obs.2,
         d3.obs = mod.b.obs.1 + mod.b.obs.3) %>% 
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