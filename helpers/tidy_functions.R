# COLNAMES: term, b, se, z, p, ci.lb, ci.ub
# Function for extracting selection model results as tidy table
tidy_psm <- function(x) {
  # get coefficients
  b <- x[[2]]$par
  # get SEs as the sqrt of the inverse of the hessian
  se.psm <- diag(sqrt(solve(x[[2]]$hessian)))
  # make row names
  # Get moderator names
  mods <- colnames(x$XX)[-1]
  # Get selection steps
  steps <- NULL
  for (i in 2:length(x$steps)) {
    steps <- c(steps, paste0(x$steps[i-1], "<p<", x$steps[i]))
  }
  # make term column
  # add the other terms
  term <- c("variance", "intrcpt", mods, steps)
  # if it's a fixed-effects model, forget "variance"
  if(x$fe == T) term <- term[-1]
  # make a data frame
  out <- data.frame(term, b = b, se = se.psm) %>% 
    # add z-scores, p-vals, CIs
    mutate(z = b/se,
           p = 2*pnorm(abs(z), lower.tail = F),
           ci.lb = b - 1.96*se,
           ci.ub = b + 1.96*se)
  return(out)
}

# Function for extracting rma results as tidy table
tidy_rma <- function(x) {
  # get model results
  out <- with(x, data.frame(b = b, #beta = beta, 
                            se, z = zval, p = pval,
                            ci.lb, ci.ub))
  # add term names & drop row names
  out <- data.frame(term = row.names(out), out)
  row.names(out) <- NULL
  # output
  return(out)
}

# Function for extracting puniform results as tidy table
tidy_punif <- function(x) {
  # get model results
  out <- with(x, data.frame(term = "intrcpt",
                            b = est.fe,
                            se = se.fe,
                            z = zval.fe,
                            ci.lb = ci.lb.fe,
                            ci.ub = ci.ub.fe))
  out <- mutate(out, p = 2*pnorm(abs(z), lower.tail = F))
  # output
  return(out)
}

# Function for extracting robu object as tidy table
tidy_robu <- function(x) {
  out <- x$reg_table
  out <- select(out, term = labels, b = b.r, se = SE, t, dfs, p = prob,
                ci.lb = CI.L, ci.ub = CI.U)
  return(out)
}

# Scratchpad for fixing up tidy_psm with term column

