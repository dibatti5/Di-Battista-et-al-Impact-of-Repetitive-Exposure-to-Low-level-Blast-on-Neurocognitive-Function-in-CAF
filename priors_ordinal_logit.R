## Code for prior simulations informing the ordinal logit modelling of the 
## RPQ 13

# Libraries
library(rstan)
library(rstanarm)
library(rethinking)  
library(StanHeaders)

## Set up R environment
options(mc.cores = parallel::detectCores()) # Parallelize chains

par(family="serif", las=1, bty="l",
    cex.axis=1, cex.lab=1, cex.main=1,
    xaxs="i", yaxs="i", mar = c(5, 5, 3, 1))

# Betancourt diagnostics
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)

## Create a prior simulation to evaluate the range of plausible estimates in 
## the relevant outcome space

# Create stan data list
prior_dat = list(
  n = 150,
  G = 3,  # Number of groups in the data
  K = 6,  # Number of categories in the RPQ13
  outcome = rep(0, 150), # empty data frame
  group = rep(c(1,2,3), 50) #random group placement
)

# Check
str(prior_dat)


# Compile model
ord_logit_prior <- stan_model("ordinal_group_fixed_prior.stan")

# Sample from model
prior_m <- sampling(
  ord_logit_prior,
  data = prior_dat,
  iter = 2000,
  chains = 4,
  seed = 123
)

# Check model
diagnostics <- util$extract_hmc_diagnostics(prior_m)
util$check_all_hmc_diagnostics(diagnostics)

samples1 <- util$extract_expectand_vals(prior_m)
base_samples <- util$filter_expectands(samples1,
                                       c('p'),
                                       check_arrays=TRUE)
util$check_all_expectand_diagnostics(base_samples)

# Plotting function for densities of probabilities
plot_posterior_probs <- function(posterior_p, scale = 100, fill = TRUE) {
  p_scaled <- posterior_p * scale
  K <- ncol(p_scaled)
  
  par(mfrow = c(1, 1), mar = c(5, 6.5, 4, 2), las = 1,
      cex.axis = 0.9, cex.lab = 1.1, mgp = c(2.5, 1, 0))
  
  plot(NULL, xlim = c(0, scale), ylim = c(0.5, K + 1.2),
       xlab = paste0("Probability (", ifelse(scale == 100, "%", ""), ")"),
       ylab = "", yaxt = "n", xaxt = "s", main = "", bty = "n",
       xaxs = "i", yaxs = "i")
  
  title("Posterior of Category Probabilities (p)", line = 2.5, cex.main = 1.3)
  axis(2, at = 1:K, labels = paste("Category", 1:K), las = 1)
  
  for (k in 1:K) {
    d <- density(p_scaled[, k], adjust = 1.1, from = 0, to = scale)
    y_scaled <- d$y / max(d$y) * 0.8 + k
    
    if (fill) {
      polygon(c(d$x, rev(d$x)), c(y_scaled, rep(k, length(d$x))),
              col = rgb(0.3, 0.5, 0.9, 0.2), border = NA)
    }
    
    lines(d$x, y_scaled, lwd = 2)
    lines(range(d$x), rep(k, 2), col = "gray60")
  }
}

# Plot the parameter estimates of the probabilities
posterior <- extract(prior_m)  # rstan 

plot_posterior_probs(posterior$p, scale = 100, fill = TRUE)

# Plot the PPC probabilities
get_ppc_probs <- function(y_pred, K) {
  n_iter <- nrow(y_pred)
  ppc <- matrix(NA, nrow = n_iter, ncol = K)
  
  for (i in 1:n_iter) {
    counts <- table(factor(y_pred[i, ], levels = 1:K))
    ppc[i, ] <- as.numeric(counts) / ncol(y_pred)
  }
  
  return(ppc)  # returns [iterations x K] matrix
}

ppc_probs <- get_ppc_probs(posterior$y_pred, 6)
plot_posterior_probs(ppc_probs)

# This produces reasonable priors assuming a relatively even split on 
# responses.

