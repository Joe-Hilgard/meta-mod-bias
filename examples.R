# examples

# let's use dataMA() to simulate a meta-analytic dataset
dat <- dataMA(k = 50,
              QRP = 0,
              sel =0,
              sigma = 0,
              min = 50, max = 500)

# let's look at it
head(dat)

# let's fit our first meta-analysis
ma1 <- rma(yi = d, sei = se, data = dat)
summary(ma1)
funnel(ma1)
trimfill(ma1)
funnel(trimfill(ma1))

# simulate null effect dataset
# let's use dataMA() to simulate a meta-analytic dataset
dat.null <- dataMA(k = 50,
              QRP = 0,
              sel =0,
              sigma = 0,
              min = 50, max = 500,
              meanD = 0)
ma2 = rma(yi = d, sei = se, data = dat.null)
summary(ma2)
funnel(ma2)
forest(ma2)

ma3 = rma(yi = d, sei = se, data = head(dat.null))
summary(ma3)
forest(ma3)

dat.bias = dataMA(k = 50,
                  QRP = 0,
                  sel =1, propB = .9,
                  sigma = 0,
                  min = 50, max = 500,
                  meanD = 0)
ma4 = rma(yi = d, sei = se, data = dat.bias)
summary(ma4)
funnel(ma4)
funnel(ma4, refline = 0)
trimfill(ma4)
funnel(trimfill(ma4), refline = 0)


# moderator effect (meta-regression)
dat.mod <- data.frame(dat, "mod" = runif(50))
ma5 = rma(yi = d, sei = se, data = dat.mod,
          mods = mod)
summary(ma5)
plot(ma5)

ma6 = rma(yi = d, sei = se, data = dat.bias, mods = se)
plot(ma6)
funnel(ma4)

#arguments for funnel() to get shaded bands for p < .10 and p < .05.
#also possible to overlay PET-PEESE regression lines not worth it at the moment. 

funnel(level=c(90, 95, 99), shade = c("white", "grey75", "grey60"), refline=0)