### Script to model mental health symptomology using Ordinal Logit Models ###

# Libraries
library(rstan)
library(rstanarm)
library(rethinking)  
library(StanHeaders)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gt)
library(gtsummary)

## Set up R environment
options(mc.cores = parallel::detectCores()) # Parallelize chains

par(family="serif", las=1, bty="l",
    cex.axis=1, cex.lab=1, cex.main=1,
    xaxs="i", yaxs="i", mar = c(5, 5, 3, 1))

# Betancourt diagnostics
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)

# Load file (Not available for public release)
df <- read.csv('data_df.csv', stringsAsFactors = FALSE)
colnames(df)

## Change mental health scores to 3-level integers ##

# PTSD
ptsd_score <- df$PTSD_TotalSeverityScore_FinalAnalysis

# Initialize as 1 for zero values
df$PTSD_cat <- ifelse(ptsd_score == 0, 1, NA)

# Get non-zero values
non_zero_ptsd <- ptsd_score[ptsd_score != 0]

# Calculate the median of non-zero values for splitting
median_non_zero_ptsd <- median(non_zero_ptsd)

# Recode non-zero values: below/equal median → 2, above → 3
df$PTSD_cat[ptsd_score > 0 & ptsd_score <= median_non_zero_ptsd] <- 2
df$PTSD_cat[ptsd_score > median_non_zero_ptsd] <- 3

# BSI Soma
bsis_score <- df$BSI_SOMA

# Initialize as 1 for zero values
df$bsis_cat <- ifelse(bsis_score == 0, 1, NA)

# Get non-zero values
non_zero_bsis <- bsis_score[bsis_score != 0]

# Calculate the median of non-zero values for splitting
median_non_zero_bsis <- median(non_zero_bsis)

# Recode non-zero values: below/equal median → 2, above → 3
df$bsis_cat[bsis_score > 0 & bsis_score <= median_non_zero_bsis] <- 2
df$bsis_cat[bsis_score > median_non_zero_bsis] <- 3

# BSI Anxiety
bsia_score <- df$BSI_ANX

# Initialize as 1 for zero values
df$bsia_cat <- ifelse(bsia_score == 0, 1, NA)

# Get non-zero values
non_zero_bsia <- bsia_score[bsia_score != 0]

# Calculate the median of non-zero values for splitting
median_non_zero_bsia <- median(non_zero_bsia)

# Recode non-zero values: below/equal median → 2, above → 3
df$bsia_cat[bsia_score > 0 & bsia_score <= median_non_zero_bsia] <- 2
df$bsia_cat[bsia_score > median_non_zero_bsia] <- 3

# BSI Depression
bsid_score <- df$BSI_DEPR

# Initialize as 1 for zero values
df$bsid_cat <- ifelse(bsid_score == 0, 1, NA)

# Get non-zero values
non_zero_bsid <- bsid_score[bsid_score != 0]

# Calculate the median of non-zero values for splitting
median_non_zero_bsid <- median(non_zero_bsid)

# Recode non-zero values: below/equal median → 2, above → 3
df$bsid_cat[bsid_score > 0 & bsid_score <= median_non_zero_bsid] <- 2
df$bsid_cat[bsid_score > median_non_zero_bsid] <- 3

## Modelling ##

# Prior modelling #

# Create stan data list for prior model
dat_prior <- list(
  G = length(unique(df$group)),
  outcome = as.integer(df$PTSD_cat) * 0, 
  group = as.integer(factor(df$group)),
  n = nrow(df),
  K = 3)

# Compile prior model
ord_logit_prior <- stan_model("ordinal_group_fixed_prior.stan")

# Sample from the prior
m_prior <- sampling(ord_logit_prior, data = dat_prior, 
                     iter = 2000, chains = 4, seed = 123)

# Real Data modelling 

# Make military controls into the first group, so they will be fixed as zero
# in the ordinal model
df$group <- factor(df$group,levels = c("Military Control","Breacher", "Sniper"))
levels(df$group)

# Create stan data lists

# PTSD
dat_ptsd <- list(
  G = length(unique(df$group)),
  outcome = as.integer(df$PTSD_cat), 
  group = as.integer(factor(df$group)),
  n = nrow(df),
  K = 3)

# Check 
str(dat_ptsd)

# BSI-S
dat_bsis <- list(
  G = length(unique(df$group)),
  outcome = as.integer(df$bsis_cat), 
  group = as.integer(factor(df$group)),
  n = nrow(df),
  K = 3)

# Check 
str(dat_bsis)

# BSI-A
dat_bsia <- list(
  G = length(unique(df$group)),
  outcome = as.integer(df$bsia_cat), 
  group = as.integer(factor(df$group)),
  n = nrow(df),
  K = 3)

# Check 
str(dat_bsia)

# BSI-D
dat_bsid <- list(
  G = length(unique(df$group)),
  outcome = as.integer(df$bsid_cat), 
  group = as.integer(factor(df$group)),
  n = nrow(df),
  K = 3)

# Check 
str(dat_bsid)

# List each stan dataset
dat_mh <- list(dat_ptsd, dat_bsis, dat_bsia, dat_bsid)

# Compile model
ord_logit <- stan_model("ordinal_group_fixed.stan")

# Sample from data model
m_mh<- list()
for (i in 1:4){
  m_mh[[i]] <- sampling(ord_logit, data = dat_mh[[i]], 
                         iter = 2000, 
                         chains = 4)}

# Check
m_mh[[1]]

# Check model diagnostics
diagnostics <- list()
check <- list()
for (i in 1:length(m_mh)){
  diagnostics[[i]] <- util$extract_hmc_diagnostics(m_mh[[i]])
  check[[i]] <- util$check_all_hmc_diagnostics(diagnostics[[i]])}

## Extract samples across the list and evaluate group differences on cumulative
## probabilities. This only works for a 2 group situation.
m_mh_posterior <- lapply(m_mh, function(x) rstan::extract(x))

# Name the lists by the mental 
names(m_mh) <- c("PCL5", "BSI-S", "BSI-A", "BSI-D")

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

# Run function across list
cdf_list <- lapply(m_mh_posterior, function(x) {
  get_param_probs_by_group(x$gamma, x$cut_points)
})
# Name the list elements by questionnaire
names(cdf_list) <- c("PCL5", "BSI-S", "BSI-A", "BSI-D")

# Create group labels to label sublist elements
group_labels <- c("Military Control", "Breacher", "Sniper")

# Create data frames and add question and group info
for (i in seq_along(cdf_list)) {
  question_name <- names(cdf_list)[i]
  
  # Convert each group's matrix to data frame and label
  cdf_list[[i]] <- lapply(seq_along(cdf_list[[i]]), function(j) {
    df <- as.data.frame(cdf_list[[i]][[j]])
    df$questionnaire <- question_name
    df$group <- group_labels[j]
    return(df)
  })
}

# Take question level dataframes and combine into a single dataframe
cdf_df <- do.call(rbind, unlist(cdf_list, recursive = FALSE))

# Change colnames
colnames(cdf_df)[1:3] <- c("None", "Low", "Moderate to High")

# Change to probs
cdf_df[c(1:3)] <- cdf_df[c(1:3)] * 100 

# Prep priors for plot
post_m_prior <- rstan::extract(m_prior)
prior_list <- get_param_probs_by_group(
  post_m_prior$gamma, post_m_prior$cut_points)

# mMake into dataframe
prior_df <- data.frame(do.call(rbind, prior_list))
prior_df[c(1:3)] <- prior_df[c(1:3)] * 100  # Convert to percentages

# Create group
prior_df$group <- rep(group_labels, each = nrow(prior_df) / length(group_labels))

# Change colnames
colnames(prior_df)[1:3] <- c("None", "Low", "Moderate to High")

### Plots ###

# Reorder factor levels for consistent coloring
cdf_df$group <- factor(cdf_df$group, levels = c("Military Control", "Breacher", "Sniper"))

# Save PDF
pdf('fig4.pdf', width = 10, height = 12, pointsize = 12)

# Layout parameters
par(
  mfrow = c(4, 3),           # 4 rows × 3 columns = 12 panels
  mar = c(3, 3, 2, 1),       # panel margins: bottom, left, top, right
  mgp = c(2, 0.6, 0),        # axis title, label, tick spacing
  oma = c(4, 5, 4, 1),       # outer margins: bottom, left, top, right
  family = "sans", 
  las = 1, bty = "l"
)

# Prep
questions <- unique(cdf_df$questionnaire)
category_labels <- c("No symptoms", "Low symptoms", "Moderate to High symptoms")
group_names <- c("Military Control", "Breacher", "Sniper")
group_colors <- c("forestgreen", "darkorange", "blue3")

# Panel loop
panel_count <- 0
for (q_idx in seq_along(questions)) {
  question <- questions[q_idx]
  post_q <- subset(cdf_df, questionnaire == question)
  
  for (k in 1:3) {
    panel_count <- panel_count + 1
    
    cat_col <- colnames(cdf_df)[k]
    post_draws <- lapply(group_names, function(gr) subset(post_q, group == gr)[[cat_col]])
    
    all_vals <- unlist(post_draws)
    xlim <- range(all_vals)
    ylim <- c(0, max(sapply(post_draws, function(x) max(density(x)$y))) * 1.2)
    
    # Plot posterior densities
    plot(density(post_draws[[1]]), lwd = 3, col = group_colors[1],
         xlim = xlim, ylim = ylim, main = "", xlab = "", ylab = "")
    lines(density(post_draws[[2]]), lwd = 3, col = group_colors[2])
    lines(density(post_draws[[3]]), lwd = 3, col = group_colors[3])
    
    # Add panel letter (A–L)
    panel_letter <- LETTERS[panel_count]
    text(x = par("usr")[1] + 0.08 * diff(par("usr")[1:2]),
         y = par("usr")[4] - 0.025 * diff(par("usr")[3:4]),
         labels = panel_letter, font = 2, cex = 1.2)
    
    # Add legend only in final panel (L)
    if (panel_count == 12) {
      legend("topright", box.lty = 0, cex = 1.2,
             legend = c("Military Control", "Breacher", "Sniper"),
             col = group_colors, lwd = 3, lty = 1)
    }
    
    # Add centered row label above middle column of each row
    if (k == 2) {
      mtext(question, side = 3, line = 1, cex = 0.8, font = 2)
    }
  }
}

# Global axis labels
mtext("Posterior Probability (%)", side = 1, font = 2, outer = TRUE, line = 2.5)
mtext("Posterior Density", side = 2, font = 2, outer = TRUE, line = 2.8, las = 0)
# Add symptom severity labels centered over each column
mtext("No symptoms", side = 3, outer = TRUE, at = 0.17, line = 2, font = 2, cex = 1.1)
mtext("Low symptoms", side = 3, outer = TRUE, at = 0.5, line = 2, font = 2, cex = 1.1)
mtext("Moderate to High symptoms", side = 3, outer = TRUE, at = 0.83, font = 2, line = 2, cex = 1.1)

# Close PDF
dev.off()

### PPCs ###

# Open PDF
pdf('fig4_ppc.pdf', width = 10, height = 12, pointsize = 12)

# Plot layout
par(
  mfrow = c(4, 3),
  mar = c(3, 3, 2, 1),
  mgp = c(2, 0.6, 0),
  oma = c(4, 4, 2, 1),
  family = "serif",
  las = 1, bty = "l"
)

questions <- c("PCL5", "BSI-S", "BSI-A", "BSI-D")
category_labels <- c("No symptoms", "Low symptoms", "Moderate to High symptoms")
group_names <- c("Military Control", "Breacher", "Sniper")
group_colors <- c("forestgreen", "darkorange", "blue3")

# Create list of PPC matrices (one per questionnaire)
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
ppc_list <- list(
  PCL5 = get_ppc_probs_by_group(m_mh_posterior[[1]]$y_pred, as.integer(df$group), 3, 3),
  `BSI-S` = get_ppc_probs_by_group(m_mh_posterior[[2]]$y_pred, as.integer(df$group), 3, 3),
  `BSI-A` = get_ppc_probs_by_group(m_mh_posterior[[3]]$y_pred, as.integer(df$group), 3, 3),
  `BSI-D` = get_ppc_probs_by_group(m_mh_posterior[[4]]$y_pred, as.integer(df$group), 3, 3)
)

# Get observed proportions
get_obs_props <- function(data, outcome, G = 3, K = 3) {
  group_ids <- as.integer(factor(data$group, levels = group_names))
  obs_mat <- matrix(0, nrow = G, ncol = K)
  for (g in 1:G) {
    y_g <- data[group_ids == g, outcome]
    counts <- table(factor(y_g, levels = 1:K))
    obs_mat[g, ] <- counts / sum(counts)
  }
  return(obs_mat)
}

obs_list <- list(
  PCL5 = get_obs_props(df, "PTSD_cat"),
  `BSI-S` = get_obs_props(df, "bsis_cat"),
  `BSI-A` = get_obs_props(df, "bsia_cat"),
  `BSI-D` = get_obs_props(df, "bsid_cat")
)

# Loop through questions and categories
panel_count <- 0
for (q in questions) {
  ppc_q <- ppc_list[[q]]
  obs_q <- obs_list[[q]]
  
  for (k in 1:3) {
    panel_count <- panel_count + 1
    
    # Extract draws for category k across groups
    ppc_k <- lapply(ppc_q, function(mat) mat[, k] * 100)
    
    # Get limits
    all_vals <- unlist(ppc_k)
    xlim <- range(all_vals)
    ylim <- c(0, max(sapply(ppc_k, function(x) max(density(x)$y))) * 1.2)
    
    # Plot first group
    plot(density(ppc_k[[1]]), lwd = 3, col = group_colors[1],
         xlim = xlim, ylim = ylim, xlab = "", ylab = "", main = "")
    lines(density(ppc_k[[2]]), lwd = 3, col = group_colors[2])
    lines(density(ppc_k[[3]]), lwd = 3, col = group_colors[3])
    
    # Add vertical lines for observed proportions
    for (g in 1:3) {
      abline(v = obs_q[g, k] * 100, lty = 4, col = group_colors[g])
    }
    
    # Axis labels
    if (k == 1) {
      mtext(side = 2, line = 3.5, at = mean(ylim), q, cex = 0.9)
    }
    if (q == questions[1]) {
      mtext(side = 3, line = 1, category_labels[k], cex = 1)
    }
    
    # Add legend only to bottom-right panel
    if (panel_count == 12) {
      legend("topright", box.lty = 0, cex = 1.0,
             legend = group_names, col = group_colors, lwd = 3, lty = 1)
    }
  }
}

# Global axis label
mtext("Posterior Probability (%)", side = 1, outer = TRUE, line = 2.5)

# Close device
dev.off()


## Table creation to facilitate results write up ##
# Convert symptom probabilities to percentages #

cdf_clean <- cdf_df %>%
  rename(
    `No symptoms` = `None`,
    `Low symptoms` = `Low`,
    `Moderate to High symptoms` = `Moderate to High`
  ) 
# Add draw ID
cdf_clean <- cdf_clean %>%
  group_by(group, questionnaire) %>%
  mutate(draw = row_number()) %>%
  ungroup()

# Pivot wider to format each questionnaire as a separate set of columns
cdf_wide <- cdf_clean %>%
  pivot_wider(
    id_cols = c(draw, group),
    names_from = questionnaire,
    values_from = `No symptoms`:`Moderate to High symptoms`,
    names_sep = " - "
  ) %>%
  relocate(group)

# Convert group to factor for order
cdf_wide$group <- factor(cdf_wide$group, levels = c("Military Control", "Breacher", "Sniper"))

# Relabel according to cut offs
colnames(cdf_wide)[-c(1,2)] <- c("PCL5 = 0", "BSI-S = 0", "BSI-A = 0", "BSI-D = 0",
                                 "0 < PCL5 ≤ 5", "0 < BSI-S ≤ 2", "0 < BSI-A ≤ 3.5",
                                 "0 < BSI-D ≤ 2", "PCL > 5", "BSI-S > 2", 
                                 "BSI-A > 3.5", "BSI-D > 2")

# Create summary table
main_tbl <- tbl_summary(
  cdf_wide %>% select(-draw),
  by = group,
  missing = "no",
  digits = all_continuous() ~ 1,
  statistic = all_continuous() ~ "{mean} ({p5} - {p95})"
) %>%
  modify_header(label = "**Symptom - Questionnaire**", all_stat_cols() ~ "**{level}**") %>%
  modify_table_styling(columns = everything(), footnote = NA) %>%
  bold_labels()

# Define contrast function
pprob_func <- function(x) mean(x > 0) * 100

# Create group-specific data
mc <- cdf_wide %>% filter(group == "Military Control") %>% select(-group, -draw)
br <- cdf_wide %>% filter(group == "Breacher") %>% select(-group, -draw)
sn <- cdf_wide %>% filter(group == "Sniper") %>% select(-group, -draw)

# Compute contrasts
bvm <- mc - br
svm <- mc - sn
bvs <- br - sn

# Summary tables for contrasts
bvm_tbl <- tbl_summary(bvm, statistic = all_continuous() ~ "{mean} ({p5} - {p95}; pprob = {pprob_func})") %>%
  modify_header(label = "**Symptom - Questionnaire**", stat_0 = "**MC – Breacher**") %>%
  modify_table_styling(columns = everything(), footnote = NA)

svm_tbl <- tbl_summary(svm, statistic = all_continuous() ~ "{mean} ({p5} - {p95}; pprob = {pprob_func})") %>%
  modify_header(label = "**Symptom - Questionnaire**", stat_0 = "**MC – Sniper**") %>%
  modify_table_styling(columns = everything(), footnote = NA)

bvs_tbl <- tbl_summary(bvs, statistic = all_continuous() ~ "{mean} ({p5} - {p95}; pprob = {pprob_func})") %>%
  modify_header(label = "**Symptom - Questionnaire**", stat_0 = "**Breacher – Sniper**") %>%
  modify_table_styling(columns = everything(), footnote = NA)

# Combine into one contrast table
contrast_tbl <- tbl_merge(list(bvm_tbl, svm_tbl, bvs_tbl))

final_tbl <- tbl_merge(
  list(main_tbl, contrast_tbl),
  tab_spanner = c("**Estimated Symptom Probabilities (%)**", "**Differences Between Groups (%)**")
) %>%
  as_gt() %>%
  gt::tab_row_group(
    group = "No symptoms",
    rows = 1:4
  ) %>%
  gt::tab_row_group(
    group = "Low symptoms",
    rows = 5:8
  ) %>%
  gt::tab_row_group(
    group = "Moderate to High symptoms",
    rows = 9:12
  ) %>%
  tab_source_note(md("The table displays the estimated percentage of participants in each symptom category by group, and the posterior contrasts between groups.")) %>%
  tab_source_note(md("Values represent posterior means and 90% credible intervals. Posterior probabilities reflect the proportion of draws above 0.")) %>%
  tab_source_note(md("Summary statistics derived from 2000 posterior draws."))

# Save as HTML
gtsave(final_tbl, "sup_tbl1.html")




