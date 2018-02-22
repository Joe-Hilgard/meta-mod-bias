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
inspectMA <- function(dataset) {
  # basic model
  rmamod <- rma(yi = d, sei = se, data = dataset)
  # test for moderator
  moderation.test <- rma(yi = d, sei = se, 
                         mods = ~id, data = dataset)
  # test for small-study effect
  egger.PET.test <- rma(yi = d, sei = se,
                        mods = ~se, data = dataset)
  # test for moderator after adjustment for small-study
  joint.test.additive <- rma(yi = d, sei = se,
                             mods = ~id + se, data = dataset)
  joint.test.interactive <- rma(yi = d, sei = se,
                                mods = ~id * se, data = dataset)
  # not sure how to interpret these parameters in the interactive model
  # test for moderator in hedges & vevea weight model
  # return test results
  out = data.frame(d.obs = summary(rmamod)$b[1],
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
                   # Egger / PET
                   d.obs.pet = summary(egger.PET.test)$b[1],
                   p.pet = summary(egger.PET.test)$pval[1],
                   b.egger = summary(egger.PET.test)$b[2],
                   p.egger = summary(egger.PET.test)$pval[2],
                   # joint PET-RMA tests
                   #Additive model 
                   joint.add.b1 = summary(joint.test.additive)$b[1],
                   joint.add.b2 = summary(joint.test.additive)$b[2],
                   joint.add.b3 = summary(joint.test.additive)$b[3],
                   joint.add.b.egger = summary(joint.test.additive)$b[4],
                   joint.add.p1 = summary(joint.test.additive)$pval[1],
                   joint.add.p2 = summary(joint.test.additive)$pval[2],
                   joint.add.p3 = summary(joint.test.additive)$pval[3],
                   joint.add.p.egger = summary(joint.test.additive)$pval[4],
                   #Interactive model
                   joint.inter.b1 = summary(joint.test.interactive)$b[1],
                   joint.inter.b2 = summary(joint.test.interactive)$b[2],
                   joint.inter.b3 = summary(joint.test.interactive)$b[3],
                   joint.inter.b1.egger = summary(joint.test.interactive)$b[4],
                   joint.inter.b2.egger = summary(joint.test.interactive)$b[5],
                   joint.inter.b3.egger = summary(joint.test.interactive)$b[6],
                   joint.inter.p1 = summary(joint.test.interactive)$pval[1],
                   joint.inter.p2 = summary(joint.test.interactive)$pval[2],
                   joint.inter.p3 = summary(joint.test.interactive)$pval[3],
                   joint.inter.p1.egger = summary(joint.test.interactive)$pval[4],
                   joint.inter.p2.egger = summary(joint.test.interactive)$pval[5],
                   joint.inter.p3.egger = summary(joint.test.interactive)$pval[6]
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
    summarize_at(.vars = vars(d.obs, mod.b.obs.1:mod.b.obs.3, d.obs.pet, 
                              joint.add.b1:joint.add.b3, joint.inter.b1:joint.inter.b3),
                 .funs = funs(mean)) %>% 
    # make cell means (NOTE: ASSUMES DUMMY CODING)
    # Might be a simpler way... lsmeans package?
    # use transmute to drop all other columns
    transmute(d.obs,
              d.obs.pet,
              d1.obs = mod.b.obs.1,
              d2.obs = mod.b.obs.1 + mod.b.obs.2,
              d3.obs = mod.b.obs.1 + mod.b.obs.3, 
              joint.add.d1.obs = joint.add.b1,
              joint.add.d2.obs = joint.add.b1 + joint.add.b2,
              joint.add.d3.obs = joint.add.b1 + joint.add.b3,
              joint.inter.d1.obs = joint.inter.b1,
              joint.inter.d2.obs = joint.inter.b1 + joint.inter.b2,
              joint.inter.d3.obs = joint.inter.b1 + joint.inter.b3)

  # Get power (or Type I) rates for each p-value test
  x.pow <- x %>% 
    summarize_at(.vars = vars(d.p, mod.p.1:mod.p.3, 
                              joint.add.p1:joint.add.p3,
                              joint.inter.p1:joint.inter.p3),
                 .funs = funs(mean(is_sig(.))))
  
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
    mutate(d1.obs = mod.b.obs.1,
           d2.obs = mod.b.obs.1 + mod.b.obs.2,
           d3.obs = mod.b.obs.1 + mod.b.obs.3) %>% 
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
  bind_rows(output.nobias, output.med_pubbias) %>% 
    gather(key, value, mod.b.obs.2:mod.b.obs.3) %>% 
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
