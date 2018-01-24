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

# original way
cool(bake(knead(mix(dough, "lightly"), "as little as possible"), temp = 350))
# piping way
dough %>% 
  mix("lightly") %>% 
  knead("as little as possible") %>% 
  bake(temp = 350) %>% 
  cool()

mtcars %>% 
  select(-gear) %>% 
  filter(mpg < 20) %>% 
  mutate(automatic = ifelse(am == 0, "manual", "automatic")) %>% 
  lm(mpg ~ wt, data = .)

# moderator effect (meta-regression)
dat.mod <- data.frame(dat, "mod" = runif(50))
ma5 = rma(yi = d, sei = se, data = dat.mod,
          mods = mod)
summary(ma5)
plot(ma5)

ma6 = rma(yi = d, sei = se, data = dat.bias, mods = se)
plot(ma6)
funnel(ma4)
