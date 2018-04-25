# Analyse results

# Load packages and functions
library(MASS)
library(tidyverse)
library(truncnorm)
library(truncdist)
library(pwr)
library(compiler)
library(metafor)
library(forcats)
library(weightr)
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
summarize_run(res.smallfx.medPB) # 
# observed d = .14, .33, .47; Type I = 83%; power = 89%, 100%; 2.2%, 35%; egger power ?%

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

summarize_run(res.hiPB)
summarize_run(res.smallfx.hiPB)
summarize_run(res.medPB.hiQRP)
summarize_run(res.smallfx.nobias)
summarize_run(res.smallfx.medPB.hiQRP)
summarize_run(res.smallfx.hiPB)

# Next steps:
# increase nSim
# inspect intermediate steps for quality, e.g. unit tests
# explore different values of vector d
# check distribution of N -- is it appropriate?
#hist(meta1$N)
# hist(meta1$N)

# implement PET-RMA model
delta <- c(0, .3, .6)
mySummary <- function(x) {
  x %>% 
    summarize(
      # Mean error
      me.b1 = mean(mod.obs.b1 - delta[1]),
      me.b2 = mean(mod.obs.b2 - delta[2]),
      me.b3 = mean(mod.obs.b3 - delta[3]),
      me.b1.add = mean(joint.add.b1 - delta[1]),
      me.b2.add = mean(joint.add.b2 - delta[2]),
      me.b3.add = mean(joint.add.b3 - delta[3]),
      me.b1.inter = mean(joint.inter.b1 - delta[1]),
      me.b2.inter = mean(joint.inter.b2 - delta[2]),
      me.b3.inter = mean(joint.inter.b3 - delta[3]),
      me.b1.sel = mean(sel.b1 - delta[1]),
      me.b1.sel = mean(sel.b2 - delta[2]),
      me.b1.sel = mean(sel.b3 - delta[3]),
      # RMSE        
      rmse.b1 = sqrt(mean((mod.obs.b1 - delta[1])^2)),
      rmse.b2 = sqrt(mean((mod.obs.b2 - delta[2])^2)),
      rmse.b3 = sqrt(mean((mod.obs.b3 - delta[3])^2)),
      rmse.b1.add = sqrt(mean((joint.add.b1 - delta[1])^2)),
      rmse.b2.add = sqrt(mean((joint.add.b2 - delta[2])^2)),
      rmse.b3.add = sqrt(mean((joint.add.b3 - delta[3])^2)),
      rmse.b1.inter = sqrt(mean((joint.inter.b1 - delta[1])^2)),
      rmse.b2.inter = sqrt(mean((joint.inter.b2 - delta[2])^2)),
      rmse.b3.inter = sqrt(mean((joint.inter.b3 - delta[3])^2)),
      rmse.b1.sel = sqrt(mean((sel.b1 - delta[1])^2)),
      rmse.b2.sel = sqrt(mean((sel.b2 - delta[2])^2)),
      rmse.b3.sel = sqrt(mean((sel.b3 - delta[3])^2)),
      # Power / Type I  
      pow.b1 = mean(mod.p1 < .05),
      pow.b2 = mean(mod.p2 < .05),
      pow.b3 = mean(mod.p3 < .05),
      pow.b1.add = mean(joint.add.p1 < .05),
      pow.b2.add = mean(joint.add.p2 < .05),
      pow.b3.add = mean(joint.add.p3 < .05),
      pow.b1.inter = mean(joint.inter.p1 < .05),
      pow.b2.inter = mean(joint.inter.p2 < .05),
      pow.b3.inter = mean(joint.inter.p3 < .05),
      pow.b1.sel = mean(sel.p1 < .05),
      pow.b2.sel = mean(sel.p2 < .05),
      pow.b3.sel = mean(sel.p3 < .05)
    )
}

final <- bind_rows(medFX_noBias = mySummary(res.nobias),
          #medFX_medPBhiQRP = mySummary(res.medPB.hiQRP),
          smallFX_noBias = mySummary(res.smallfx.nobias),
          smallFX_medPBhiQRP = mySummary(res.smallfx.medPB.hiQRP),
          .id = "id")
write.csv(final, "final_output.csv", row.names = F)


# Plotting for poster ----
theme_poster <- theme(text = element_text(size = 32),
                      legend.position = "none")

temp <- bind_rows(noBias = res.nobias, 
                  bias = res.medPB.hiQRP, 
                  .id = "id") %>% 
  select(id,
         mod.obs.b1:mod.obs.b3, 
         joint.add.b1:joint.add.b3, 
         joint.inter.b1:joint.inter.b3,
         sel.b1:sel.b3) %>% 
  gather(key = "key", value = "value", mod.obs.b1:sel.b3) %>% 
  separate(key, into = c("junk", "estimator", "coefficient"), sep = "\\.") %>% 
  # augment with true values & fix factor order
  mutate(delta = ifelse(coefficient == "b1", 0,
                        ifelse(coefficient == "b2", 0.3, 0.6)),
         rmse = sqrt((value - delta)^2),
         id = fct_relevel(id, c("noBias", "bias")),
         estimator = fct_relevel(estimator, c("obs", "add", "inter")))

temp2 <- bind_rows(noBias = res.smallfx.nobias, 
                  bias = res.smallfx.medPB.hiQRP, 
                  .id = "id") %>% 
  select(id,
         mod.obs.b1:mod.obs.b3, 
         joint.add.b1:joint.add.b3, 
         joint.inter.b1:joint.inter.b3,
         sel.b1:sel.b3) %>% 
  gather(key = "key", value = "value", mod.obs.b1:sel.b3) %>% 
  separate(key, into = c("junk", "estimator", "coefficient"), sep = "\\.") %>% 
  # augment with true values
  mutate(delta = ifelse(coefficient == "b1", 0,
                        ifelse(coefficient == "b2", 0.2, 0.4)),
         rmse = sqrt((value - delta)^2),
         id = fct_relevel(id, c("noBias", "bias")),
         estimator = fct_relevel(estimator, c("obs", "add", "inter")))

# plot means
temp %>% 
  group_by(id, estimator, coefficient, delta) %>% 
  summarize(m = mean(value),
            q.lower = quantile(value, .025),
            q.upper = quantile(value, .975)) %>% 
  ggplot(aes(x = coefficient, y = m, shape = estimator)) +
  geom_hline(aes(yintercept = delta), col = 'grey25') +
  geom_pointrange(aes(ymin = q.lower, ymax = q.upper),
                  size = 2, position = position_dodge(width = .5)) +
  facet_grid(~id) +
  scale_y_continuous("Estimate (d)",
                     breaks = c(-0.3, 0, 0.3, 0.6, 0.9)) +
  theme_poster
ggsave("me_bigfx.png", width = 10.25, height = 4)

temp2 %>% 
  group_by(id, estimator, coefficient, delta) %>% 
  summarize(m = mean(value),
            q.lower = quantile(value, .025),
            q.upper = quantile(value, .975)) %>% 
  ggplot(aes(x = coefficient, y = m, shape = estimator)) +
  geom_hline(aes(yintercept = delta), col = 'grey25') +
  geom_pointrange(aes(ymin = q.lower, ymax = q.upper),
                  size = 2, position = position_dodge(width = .5)) +
  facet_grid(~id) +
  scale_y_continuous("Estimate (d)",
                     breaks = c(-0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8)) +
  theme_poster
ggsave("me_smfx.png", width = 10.25, height = 4)


# plot RMSE
# temp %>% 
#   group_by(id, estimator, coefficient) %>% 
#   summarize(m = mean(rmse),
#             q.lower = quantile(rmse, .025),
#             q.upper = quantile(rmse, .975)) %>% 
#   ggplot(aes(x = coefficient, y = m, shape = estimator,
#              ymin = q.lower, ymax = q.upper)) +
#   geom_pointrange(size = 2, position = position_dodge(width = .5)) +
#   facet_wrap(~id) +
#   scale_y_continuous("RMSE") +
#   theme_poster
# ggsave("rmse_bigfx.png", width = 10.25, height = 4)
# 
# temp2 %>% 
#   group_by(id, estimator, coefficient) %>% 
#   summarize(m = mean(rmse),
#             q.lower = quantile(rmse, .025),
#             q.upper = quantile(rmse, .975)) %>% 
#   ggplot(aes(x = coefficient, y = m, shape = estimator,
#              ymin = q.lower, ymax = q.upper)) +
#   geom_pointrange(size = 2, position = position_dodge(width = .5)) +
#   facet_wrap(~id) +
#   scale_y_continuous("RMSE") +
#   theme_poster
# ggsave("rmse_smfx.png", width = 10.25, height = 4)

# Power / Type I
temp.p <- bind_rows(noBias = res.nobias, 
                    bias = res.medPB.hiQRP, 
                    .id = "id") %>% 
  select(id,
         mod.p1:mod.p3, 
         joint.add.p1:joint.add.p3, 
         joint.inter.p1:joint.inter.p3,
         sel.p1:sel.p3) %>% 
  gather(key = "key", value = "value", mod.p1:sel.p3) %>% 
  separate(key, into = c("junk", "estimator", "coefficient"), sep = "\\.",
           # need fill = "left" to deal with shorter "mod.p1" vs "joint.add.p1"
           fill = "left") %>% 
  mutate(id = fct_relevel(id, "noBias", "bias"),
         estimator = fct_relevel(estimator, c("mod", "add", "inter")))

temp.p2 <- bind_rows(noBias = res.smallfx.nobias, 
                    bias = res.smallfx.medPB.hiQRP, 
                    .id = "id") %>% 
  select(id,
         mod.p1:mod.p3, 
         joint.add.p1:joint.add.p3, 
         joint.inter.p1:joint.inter.p3,
         sel.b1:sel.b3) %>% 
  gather(key = "key", value = "value", mod.p1:sel.p3) %>% 
  separate(key, into = c("junk", "estimator", "coefficient"), sep = "\\.",
           # need fill = "left" to deal with shorter "mod.p1" vs "joint.add.p1"
           fill = "left") %>% 
  mutate(id = fct_relevel(id, "noBias", "bias"),
         estimator = fct_relevel(estimator, c("mod", "add", "inter")))

temp.p %>% 
  group_by(id, estimator, coefficient) %>% 
  summarize(sig = mean(value < .05)) %>% 
  ggplot(aes(x = coefficient, y = sig, 
             shape = interaction(estimator, id),
             group = estimator)) +
  geom_point(size = 5, position = position_dodge(width = .5)) +
  scale_shape_manual(values = c(16, 17, 15, 1, 2, 0)) +
  theme_poster
ggsave("nhst_bigfx.png", width = 6, height = 4)

temp.p2 %>% 
  group_by(id, estimator, coefficient) %>% 
  summarize(sig = mean(value < .05)) %>% 
  ggplot(aes(x = coefficient, y = sig, 
             shape = interaction(estimator, id),
             group = estimator)) +
  geom_point(size = 5, position = position_dodge(width = .5)) +
  scale_shape_manual(values = c(16, 17, 15, 1, 2, 0)) +
  theme_poster
ggsave("nhst_smfx.png", width = 6, height = 4)

# temp.p2 %>% 
#   group_by(id, estimator, coefficient) %>% 
#   summarize(sig = mean(value < .05)) %>% 
#   ggplot(aes(x = coefficient, y = sig, shape = estimator)) +
#   geom_point(size = 3, position = position_dodge(width = .5)) +
#   facet_wrap(~id)





# approach 2: scatterpoints
# Looks too muddy
# temp %>% 
#   group_by(id, estimator, coefficient, delta) %>% 
#   ggplot(aes(x = coefficient, y = value, shape = estimator)) +
#   geom_boxplot(position = position_dodge(width = .5)) +
#   geom_jitter(alpha = .5,
#               position = position_jitterdodge(dodge.width = .5, jitter.width = .25)) +
#   facet_grid(~id) +
#   geom_point(aes(y = delta), col = 'darkred', size = 2) +
#   scale_y_continuous(breaks = c(-0.3, 0, 0.3, 0.6, 0.9))
