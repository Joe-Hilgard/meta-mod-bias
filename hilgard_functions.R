# This file creates four functions
# modMA() runs dataMA() three times to create a dataset with meta-moderators
# inspectMA() fits the various models to the output of modMA() & retrieves stats
# runStudy() performs modMA() and inspectMA() in a loop
# summarize_run() condenses runStudy() output into means and power rates


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
  # combine three subgroups into dataset, label subgroups with "id"
  data.out <- bind_rows(d1, d2, d3, .id = "id") %>% 
    mutate(id = factor(id))
  # output simulated dataset
  return(data.out)
}

# function to analyze results of dataMA
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
runStudy <- function(nSim, k, d, sel, propB, QRP) {
  output <- NULL; 
  for (i in 1:nSim) {
    # make a dataset of simulated studies
    # k studies each, effect sizes given by d
    tempdata <- modMA(k = k, d = d, 
                      sel = sel, propB = propB, QRP = QRP)
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
    # make cell means (NOTE: ASSUMES CONTRAST CODING)
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
