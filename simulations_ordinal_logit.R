# Simulation data for the ordered logit model

# Libraries
library(rstan)
library(rstanarm)
library(rethinking)  
library(loo)
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

## Simulate data that show group-level differences in likert responsesn across
## 3 groups

set.seed(123)

# Number of observations per group
n_per_group <- 100

# Number of groups
G <- 3

# Group-level latent locations (gamma[group])
gamma <- c(-1.0, 0.5, 1.5)

# Ordered cutpoints (must be increasing)
cutpoints <- c(-1.5, -0.5, 0, 0.5, 1.5)  # 5 cut points = 6 Likert outcomes

# Function to simulate ordered logit responses
simulate_ordered_logit <- function(eta, cutpoints) {
  # Cumulative probabilities
  probs <- numeric(length(cutpoints) + 1)
  
  probs[1] <- plogis(cutpoints[1] - eta)
  for (k in 2:length(cutpoints)) {
    probs[k] <- plogis(cutpoints[k] - eta) - sum(probs[1:(k-1)])
  }
  probs[length(probs)] <- 1 - sum(probs[1:(length(probs)-1)])
  
  # Sample category
  sample(1:length(probs), size = 1, prob = probs)
}

# Function to calculate category probabilities
ordered_logit_probs <- function(eta, cutpoints) {
  K <- length(cutpoints) + 1
  probs <- numeric(K)
  
  probs[1] <- plogis(cutpoints[1] - eta)
  for (k in 2:(K - 1)) {
    probs[k] <- plogis(cutpoints[k] - eta) - sum(probs[1:(k - 1)])
  }
  probs[K] <- 1 - sum(probs[1:(K - 1)])
  
  return(probs)
}

# Generate full dataset
group_ids <- rep(1:G, each = n_per_group)
outcomes <- mapply(function(g) simulate_ordered_logit(gamma[g], cutpoints), group_ids)

# Get probabilities for plot
group_probs <- sapply(gamma, ordered_logit_probs, cutpoints = cutpoints)
# Each column = a group, each row = a category
colnames(group_probs) <- paste("Group", 1:G)
rownames(group_probs) <- paste("Category", 1:(length(cutpoints) + 1))

# Create data frame
sim_data <- data.frame(
  group = group_ids,
  outcome = as.integer(outcomes)
)

# Check
head(sim_data)

# Visualize by group
groups <- unique(sim_data$group)
par(mfrow = c(length(groups), 1), mgp = c(3, 1.5, 0))  # one plot per row

for (g in groups) {
  util$plot_line_hist(sim_data$outcome[sim_data$group == g], 0.5, 6.5, 1,
                      xlab = paste("Ordinal - Group", g))
}

# Visualize probs by category by group
par(mfrow = c(1,1) ) 
# Transpose so each row is a group
prob_matrix <- t(group_probs)

barplot(t(prob_matrix), beside = TRUE, horiz = TRUE,
        col = gray.colors(nrow(group_probs), start = 0.9, end = 0.3),
        xlab = "Probability",
        legend.text = rownames(group_probs),
        args.legend = list(x = "topright", title = "Category"),
        main = "Theoretical Category Probabilities by Group")


# Create stan data list
sim_dat = list(
  n = nrow(sim_data),
  G = length(unique(sim_data$group)),  
  K = length(unique(sim_data$outcome)),  # Number of categories
  outcome = sim_data$outcome,
  group = sim_data$group
)

# Check
str(sim_dat)


# Compile model
ord_logit <- stan_model("ordinal_group_fixed.stan")

# Sample from model
sim_m1 <- sampling(
  ord_logit,
  data = sim_dat,
  iter = 2000,
  chains = 4,
  seed = 123
)

# Check model
diagnostics <- util$extract_hmc_diagnostics(sim_m1)
util$check_all_hmc_diagnostics(diagnostics)

samples1 <- util$extract_expectand_vals(sim_m1)
base_samples <- util$filter_expectands(samples1,
                                       c('p'),
                                       check_arrays=TRUE)
util$check_all_expectand_diagnostics(base_samples)

##  Plot probabilities ##

# Function to create the probabilities by group
get_param_probs_by_group <- function(gamma_draws, cutpoint_draws) {
  n_iter <- nrow(gamma_draws)
  G <- ncol(gamma_draws)
  K <- ncol(cutpoint_draws) + 1
  
  posterior_p_list <- vector("list", G)
  
  for (g in 1:G) {
    p_g <- matrix(NA, nrow = n_iter, ncol = K)
    
    for (s in 1:n_iter) {
      gamma_sg <- gamma_draws[s, g]
      cutpoints_s <- cutpoint_draws[s, ]
      
      cdf <- plogis(cutpoints_s - gamma_sg)
      probs <- numeric(K)
      probs[1] <- cdf[1]
      for (k in 2:(K - 1)) {
        probs[k] <- cdf[k] - cdf[k - 1]
      }
      probs[K] <- 1 - cdf[K - 1]
      
      p_g[s, ] <- probs
    }
    
    posterior_p_list[[g]] <- p_g
  }
  
  return(posterior_p_list)
}


# Plotting function for densities of probabilities
plot_group_posterior_probs <- function(posterior_p_list, scale = 100, 
                                       group_names = NULL, fill = TRUE) {
  # posterior_p_list: a list of [iterations x K] matrices, one per group
  G <- length(posterior_p_list)
  K <- ncol(posterior_p_list[[1]])
  p_scaled_list <- lapply(posterior_p_list, function(p) p * scale)
  
  # Default group names
  if (is.null(group_names)) {
    group_names <- paste("Group", seq_len(G))
  }
  
  # Colors for groups
  group_colors <- RColorBrewer::brewer.pal(G, "Set1")
  
  # Setup plotting window
  par(mfrow = c(1, 1), mar = c(5, 6.5, 4, 2), las = 1,
      cex.axis = 0.9, cex.lab = 1.1, mgp = c(2.5, 1, 0))
  
  plot(NULL, xlim = c(0, scale), ylim = c(0.5, K + 1.2),
       xlab = paste0("Probability (", ifelse(scale == 100, "%", ""), ")"),
       ylab = "", yaxt = "n", xaxt = "s", main = "", bty = "n",
       xaxs = "i", yaxs = "i")
  
  title("Posterior of Category Probabilities by Group", line = 2.5, cex.main = 1.3)
  axis(2, at = 1:K, labels = paste("Category", 1:K), las = 1)
  
  for (k in 1:K) {
    for (g in 1:G) {
      d <- density(p_scaled_list[[g]][, k], adjust = 1.1, from = 0, to = scale)
      y_scaled <- d$y / max(d$y) * 0.8 + k
      
      if (fill) {
        polygon(c(d$x, rev(d$x)), c(y_scaled, rep(k, length(d$x))),
                col = adjustcolor(group_colors[g], alpha.f = 0.2), border = NA)
      }
      lines(d$x, y_scaled, lwd = 2, col = group_colors[g])
    }
    # Horizontal baseline under each category
    lines(c(0, scale), rep(k, 2), col = "gray60")
  }
  
  legend("topright", legend = group_names, fill = adjustcolor(group_colors, alpha.f = 0.3),
         border = NA, bty = "n", title = "Group")
}

# Extract the posterior
posterior <- extract(sim_m1)  
gamma_draws <- posterior$gamma   
cutpoint_draws <- posterior$cut_points  

str(cutpoint_draws)
# Derive group-based probabilities
posterior_p_list <- get_param_probs_by_group(gamma_draws, cutpoint_draws)

# Plot posterior densities by group
plot_group_posterior_probs(posterior_p_list, group_names = c("Military Controls", "Breachers", "Snipers"))

# Plot a PPC (group agnostic) using Betancourt's code
par(mfrow=c(1, 1), mar=c(5, 5, 1, 1))

util$plot_hist_quantiles(samples1, 'y_pred', 0.5, 6.5, 1,
                         baseline_values=sim_data$outcome, xlab="Ordinal")

## PPC of probabilities by group ##
get_ppc_probs_by_group <- function(y_pred, group, G, K) {
  n_iter <- nrow(y_pred)
  posterior_p_list <- vector("list", G)
  
  for (g in 1:G) {
    # Indices of observations in group g
    group_idx <- which(group == g)
    n_g <- length(group_idx)
    
    # Initialize matrix to hold category proportions: [iterations x K]
    p_g <- matrix(NA, nrow = n_iter, ncol = K)
    
    for (s in 1:n_iter) {
      draws <- y_pred[s, group_idx]
      counts <- table(factor(draws, levels = 1:K))
      p_g[s, ] <- as.numeric(counts) / n_g
    }
    
    posterior_p_list[[g]] <- p_g
  }
  
  return(posterior_p_list)
}

ppc_probs_by_group <- get_ppc_probs_by_group(posterior$y_pred, 
                                             sim_data$group, 
                                             G = 3, K = 6)

# Derive observed proportions
get_observed_props <- function(data, group_var, outcome_var, G, K) {
  obs_props <- matrix(NA, nrow = G, ncol = K)
  for (g in 1:G) {
    counts <- table(factor(data[data[[group_var]] == g, outcome_var], levels = 1:K))
    obs_props[g, ] <- as.numeric(counts) / sum(counts)
  }
  return(obs_props)
}

obs_props <- get_observed_props(sim_data, 
                                group_var = "group", 
                                outcome_var = "outcome", 
                                G = 3, K = 6)

# Function for "facet-based" plot
plot_ppc_probs_facet <- function(ppc_probs_list, obs_props, 
                                 group_names = NULL, scale = 100) {
  G <- length(ppc_probs_list)
  K <- ncol(ppc_probs_list[[1]])
  if (is.null(group_names)) group_names <- paste("Group", 1:G)
  
  layout(matrix(1:G, nrow = G, byrow = TRUE), heights = rep(1, G))
  
  # Larger left and top margins
  par(oma = c(4, 2, 4, 2), mar = c(2, 6, 2.5, 1), las = 1,
      cex.axis = 0.95, cex.lab = 1.1, mgp = c(3, 1, 0),
      xaxs = "i", yaxs = "i")
  
  category_spacing <- 2     # More space between categories
  density_height <- 0.9     # Controls vertical stretch of each density
  max_y <- category_spacing * K + 1  # Room for top
  
  for (g in 1:G) {
    plot(NULL, xlim = c(0, scale), ylim = c(0.5, max_y),
         xlab = "", ylab = "",
         yaxt = "n", xaxt = ifelse(g == G, "s", "n"),
         main = group_names[g],
         bty = "n")
    
    axis(2,
         at = seq(category_spacing, category_spacing * K, by = category_spacing),
         labels = paste("Cat", 1:K), las = 1)
    
    if (g == G) axis(1)
    
    for (k in 1:K) {
      y_base <- category_spacing * k
      d <- density(ppc_probs_list[[g]][, k] * scale, adjust = 1.2, from = 0, to = scale)
      y_scaled <- d$y / max(d$y) * density_height + y_base
      
      polygon(c(d$x, rev(d$x)), c(y_scaled, rep(y_base, length(d$x))),
              col = rgb(0.3, 0.5, 0.9, 0.2), border = NA)
      lines(d$x, y_scaled, lwd = 2, col = "blue4")
      
      # Observed proportion
      x_obs <- obs_props[g, k] * scale
      segments(x_obs, y_base - 0.35, x_obs, y_base + 0.35, col = "black", lwd = 2)
      
      lines(c(0, scale), rep(y_base, 2), col = "gray80")
    }
  }
  
  mtext("Predicted Probability (%)", 1, outer = TRUE, line = 2.8, cex = 1.1)
  mtext("Posterior Predictive Probabilities by Group", 3, outer = TRUE, line = 1.5, cex = 1.4)
}

# Plot the PPC
plot_ppc_probs_facet(ppc_probs_by_group, 
                     obs_props, 
                     group_names = c("Military Controls", "Breachers", "Snipers"))

## Estimates are recovered well ##


