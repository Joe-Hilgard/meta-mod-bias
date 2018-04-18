# This file creates four functions
# modMA() runs simMA() three times to create a dataset with meta-moderators
# inspectMA() fits the various models to the output of modMA() & retrieves stats
# runStudy() performs modMA() and inspectMA() in a loop
# summarize_run() condenses runStudy() output into means and power rates

# There are two data-generating functions
# The default approach (used by modMA and runStudy) generates study results,
# then publishes the nonsignificant ones with some probability set by "censor"
# The alternative approach (used by modMA.post and runStudy.post) forces some
# percentage of results (selProp) to be statistically significant.
# These two approaches yield radically different degrees of publication bias.
# IMO, to get the default approach to yield realistic levels of stat. significance
# requires either strong pub bias or heavy p-hacking.

# Remember, available levels for arguments are:
# censor: "none", "low", "medium", "high"
# qrpEnv: "none", "medium", "high"

modMA <- function(k, delta, tau = 0,
                  empN = FALSE, maxN = 200, minN = 20, meanN = 50,
                  censor = "none", qrpEnv = "none", empN.boost = 0) {
  d1 <- simMA(k = k, delta = delta[1], tau = tau, 
               empN = empN, maxN = maxN, minN = minN, meanN = meanN, 
               censor = censor, qrpEnv = qrpEnv, empN.boost = empN.boost)
  d2 <- simMA(k = k, delta[2], tau = tau, 
               empN = empN, maxN = maxN, minN = minN, meanN = meanN,  
               censor = censor, qrpEnv = qrpEnv, empN.boost = empN.boost)
  d3 <- simMA(k = k, delta[3], tau = tau, 
               empN = empN, maxN = maxN, minN = minN, meanN = meanN, 
               censor = censor, qrpEnv = qrpEnv, empN.boost = empN.boost)
  # bind three datasets together into one & coerce to data.frame
  data.out <- (bind_rows(data.frame(d1), data.frame(d2), data.frame(d3), .id = "id"))
  # convert id to a factor for meta-regression
  data.out$id <- as.factor(data.out$id)
  # output simulated dataset
  return(data.out)
}

# Make a function that runs simMA three times,
# once for each population,
# then mixes them together
# TODO: consider doing something cleverer so that user can specify number of levels of "id"
modMA.post <- function(k, delta, tau = 0,
                  empN = FALSE, maxN = 200, minN = 20, meanN = 50,
                  selProp = 0, qrpEnv = "none", empN.boost = 0) {
  d1 <- dataMA(k = k, delta = delta[1], tau = tau, 
              empN = empN, maxN = maxN, minN = minN, meanN = meanN, 
              selProp = selProp, qrpEnv = qrpEnv, empN.boost = empN.boost)
  d2 <- dataMA(k = k, delta[2], tau = tau, 
              empN = empN, maxN = maxN, minN = minN, meanN = meanN,  
              selProp = selProp, qrpEnv = qrpEnv, empN.boost = empN.boost)
  d3 <- dataMA(k = k, delta[3], tau = tau, 
              empN = empN, maxN = maxN, minN = minN, meanN = meanN, 
              selProp = selProp, qrpEnv = qrpEnv, empN.boost = empN.boost)
  # bind three datasets together into one & coerce to data.frame
  data.out <- (bind_rows(data.frame(d1), data.frame(d2), data.frame(d3), .id = "id"))
  # convert id to a factor for meta-regression
  data.out$id <- as.factor(data.out$id)
  # output simulated dataset
  return(data.out)
}

# function to analyze results of simMA
# TODO: fetch SE of parameters for CI inspection?
# TODO: Do I need to figure something nicer out using map(), safely(), possibly()?
# In meta-showdown, it looks like Felix ran them as lists
#    see https://github.com/nicebread/meta-showdown/blob/master/2-analysisFramework.R#L50
# use of tryCatch: https://github.com/nicebread/meta-showdown/blob/master/MA-methods/2-p-curve.R#L19
# use of tryCatch to flip from "REML" to "DL" when necessary in PET/PEESE:
#    see https://github.com/nicebread/meta-showdown/blob/master/MA-methods/3-PET-PEESE.R#L13
# use of tryCatch to make failed weightfunct output into NULL / NA:
#    see https://github.com/nicebread/meta-showdown/blob/master/MA-methods/7-Selection%20Models.R#L12

inspectMA <- function(dataset) {
  # basic model
  rmamod <- tryCatch(
    rma(yi = d, sei = se, data = dataset),
    error = function(e) rma(yi = d, sei = se, data = dataset,
                            control = list(stepadj = .5))
    )
  # test for moderator
  moderation.test <- tryCatch(
    rma(yi = d, sei = se, mods = ~id, data = dataset),
    error = function(e) rma(yi = d, sei = se, mods = ~id, data = dataset,
                            control = list(stepadj = .5))
    )
  # test for small-study effect
  egger.PET.test <- tryCatch(
    rma(yi = d, sei = se, mods = ~se, data = dataset),
    error = function(e) rma(yi = d, sei = se, mods = ~se, data = dataset,
                            control = list(stepadj = .5))
    )
  # test for moderator after adjustment for small-study
  joint.test.additive <- rma(yi = d, sei = se,
                             mods = ~id + se, data = dataset,
                             control = list(stepadj = .5))
  joint.test.interactive <- rma(yi = d, sei = se,
                                mods = ~id * se, data = dataset,
                                control = list(stepadj = .5))
  # TODO: test for moderator in hedges & vevea weight model
  weightmodel <- with(dataset,
                      weightfunct(d, v, mods = ~id))
  # return test results
  res.basic <- summary(rmamod)
  res.mod <- summary(moderation.test)
  res.pet <- summary(egger.PET.test)
  res.add <- summary(joint.test.additive)
  res.int <- summary(joint.test.interactive)
  res.sel <- tidy_psm(weightmodel)
  
  # Collate results into data frame
  # Maybe the tryCatch way to do this would be to just tidy() all results
  #    bind them, gather them, then filter and spread them.
  out = data.frame(d.obs = res.basic$b[1],
                   se.obs = res.basic$se[1],
                   d.p = res.basic$pval[1],
                   # Moderators
                   # Note that .1 is the intercept, reference group
                   mod.obs.b1 = res.mod$b[1],
                   mod.obs.b2 = res.mod$b[2],
                   mod.obs.b3 = res.mod$b[3],
                   mod.p1 = res.mod$pval[1],
                   mod.p2 = res.mod$pval[2],
                   mod.p3 = res.mod$pval[3],
                   # Egger / PET
                   d.obs.pet = res.pet$b[1],
                   p.pet = res.pet$pval[1],
                   b.egger = res.pet$b[2],
                   p.egger = res.pet$pval[2],
                   # joint PET-RMA tests
                   #Additive model 
                   joint.add.b1 = res.add$b[1],
                   joint.add.b2 = res.add$b[2],
                   joint.add.b3 = res.add$b[3],
                   joint.add.b.egger = res.add$b[4],
                   joint.add.p1 = res.add$pval[1],
                   joint.add.p2 = res.add$pval[2],
                   joint.add.p3 = res.add$pval[3],
                   joint.add.p.egger = res.add$pval[4],
                   #Interactive model
                   joint.inter.b1 = res.int$b[1],
                   joint.inter.b2 = res.int$b[2],
                   joint.inter.b3 = res.int$b[3],
                   joint.inter.b1.egger = res.int$b[4],
                   joint.inter.b2.egger = res.int$b[5],
                   joint.inter.b3.egger = res.int$b[6],
                   joint.inter.p1 = res.int$pval[1],
                   joint.inter.p2 = res.int$pval[2],
                   joint.inter.p3 = res.int$pval[3],
                   joint.inter.p1.egger = res.int$pval[4],
                   joint.inter.p2.egger = res.int$pval[5],
                   joint.inter.p3.egger = res.int$pval[6],
                   # Hedges-Vevea model
                   sel.b1 = res.sel[2, "b"],
                   sel.b2 = res.sel[3, "b"],
                   sel.b3 = res.sel[4, "b"],
                   sel.p1 = res.sel[2, "p"],
                   sel.p2 = res.sel[3, "p"],
                   sel.p3 = res.sel[4, "p"]
  )
  return(data.frame(out))
}

# function for performing modMA() and inspectMA() nSim number of times,
# storing the output in a data frame.
runStudy <- function(nSim, k, delta, tau = 0,
                     empN = 0, maxN = 200, minN = 20, meanN = 50,
                     censor = "none", qrpEnv = "none", empN.boost = 0) {
  output <- NULL; 
  for (i in 1:nSim) {
    # make a dataset of simulated studies
    # k studies each, effect sizes given by d
    tempdata <- modMA(k, delta, tau = 0,
                      empN = empN, maxN = maxN, minN = minN, meanN = meanN,
                      censor = censor, qrpEnv = qrpEnv, empN.boost = empN.boost)
    # fit our various meta-analytic models to the simulated data
    tempresult <- inspectMA(tempdata)
    # add this iteration's results to the output object
    output <- bind_rows(output, tempresult)
  }
  output
}

# function for performing modMA.post() and inspectMA() nSim number of times,
# storing the output in a data frame.
runStudy.post <- function(nSim, k, delta, tau = 0,
                     empN = 0, maxN = 200, minN = 20, meanN = 50,
                     selProp = 0, qrpEnv = "none", empN.boost = 0) {
  output <- NULL; 
  for (i in 1:nSim) {
    # make a dataset of simulated studies
    # k studies each, effect sizes given by d
    tempdata <- modMA.post(k, delta, tau = 0,
                      empN = empN, maxN = maxN, minN = minN, meanN = meanN,
                      selProp = selProp, qrpEnv = qrpEnv, empN.boost = empN.boost)
    # fit our various meta-analytic models to the simulated data
    tempresult <- inspectMA(tempdata)
    # add this iteration's results to the output object
    output <- bind_rows(output, tempresult)
  }
  output
}

# functions for summarizing a simulation series in terms of mean estimates and
#    error rates
# TODO: expand and refine to capture ME, RMSE, power, additive and interactive models
# TODO: Might separate into two functions... a little hard on my eyes as it is with both estimates and p-vals
# TODO: is there any way to get an ANOVA-like output rather than pairwise t-test?
is_sig <- function(x) x < .05
summarize_run <- function(x) {
  # Get estimates of d, moderator parameters, bias parameter, bias-adjusted PET
  x.est <- x %>% 
    summarize_at(.vars = vars(d.obs, mod.obs.b1:mod.obs.b3, d.obs.pet, 
                              joint.add.b1:joint.add.b3, joint.inter.b1:joint.inter.b3,
                              sel.b1:sel.b3),
                 .funs = funs(mean(., na.rm = T))) %>% 
    # make cell means (NOTE: ASSUMES DUMMY CODING)
    # Might be a simpler way... lsmeans package?
    # use transmute to drop all other columns
    transmute(d.obs,
              d.obs.pet,
              # moderator model
              d1.obs = mod.obs.b1,
              d2.obs = mod.obs.b1 + mod.obs.b2,
              d3.obs = mod.obs.b1 + mod.obs.b3, 
              # PET model
              joint.add.d1.obs = joint.add.b1,
              joint.add.d2.obs = joint.add.b1 + joint.add.b2,
              joint.add.d3.obs = joint.add.b1 + joint.add.b3,
              # PET * moderator model
              joint.inter.d1.obs = joint.inter.b1,
              joint.inter.d2.obs = joint.inter.b1 + joint.inter.b2,
              joint.inter.d3.obs = joint.inter.b1 + joint.inter.b3,
              # Selection model
              sel.d1.obs = sel.b1,
              sel.d2.obs = sel.b1 + sel.b2,
              sel.d3.obs = sel.b1 + sel.b3
              )

  # Get power (or Type I) rates for each p-value test
  x.pow <- x %>% 
    summarize_at(.vars = vars(d.p, mod.p1:mod.p3, 
                              joint.add.p1:joint.add.p3,
                              joint.inter.p1:joint.inter.p3,
                              sel.p1:sel.p3),
                 .funs = funs(mean(is_sig(.), na.rm = T)))
  
  return(bind_cols(x.est, x.pow))
}

# Plot observed cell means across runs
plotCellMeans <- function(data1, data2, name1, name2) {
  # Label facets with name arguments
  data1$bias <- name1
  data2$bias <- name2
  # Generate cell means (assumes dummy coding!)
  # and plot them
  bind_rows(data1, data2) %>% 
    mutate(d1.obs = mod.obs.b1,
           d2.obs = mod.obs.b1 + mod.obs.b2,
           d3.obs = mod.obs.b1 + mod.obs.b3) %>% 
    gather(key, value, d1.obs:d3.obs) %>% 
    ggplot(aes(x = value)) +
    geom_histogram() +
    facet_grid(bias ~ key) +
    scale_x_continuous(limits = c(-.2, 1))
}

# Plot moderator's observed effect size across runs
plotParamMeans <- function(data1, data2, name1, name2) {
  # Label facets with name arguments
  data1$bias <- name1
  data2$bias <- name2
  # plot beta values
  bind_rows(data1, data2) %>% 
    gather(key, value, mod.obs.b2:mod.obs.b3) %>% 
    ggplot(aes(x = value)) +
    geom_histogram() +
    facet_grid(bias ~ key) +
    scale_x_continuous(limits = c(0, 1))  +
    ggtitle("Moderation between d = 0, 0.3, 0.6")
}

# convenient funnel in the style I like, centered at zero w/ significance bands
myFunnel <- function(x) {
  funnel(x,
         refline=0, level=c(90, 95, 99), 
         shade = c("white", "grey75", "grey60"), back = "gray90")
}

# Piping lesson
mtcars %>% 
  filter(am == 1) %>% 
  lm(mpg ~ wt, data = .)
  