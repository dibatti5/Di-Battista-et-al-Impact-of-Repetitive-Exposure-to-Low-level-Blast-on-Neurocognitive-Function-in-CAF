### Script to model Concussion symptoms using an Ordered Logit Model ###

# Libraries
library(rstan)
library(rstanarm)
library(rethinking)  
library(StanHeaders)
library(gt)
library(gtsummary)
library(ggplot2)
library(dplyr)
library(tidyr)

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

# Create a dataframe with just the RPQ 13 scores
mod_rpq <- df[c(30:45)]

# Make complete case, first change 99 to NA
mod_rpq[mod_rpq == 99] <- NA

# Change names
colnames(mod_rpq) <- c("Headache", "Dizziness", "Vomiting", 
  "Noise Sensitivity", "Sleep Disturbance",
  "Fatigue", "Irritable", "Depressed",
  "Frustrated", "Forgetful", "Poor Concentration", 
  "Taking Longer to Think", "Blurred Vision", 
  "Light Sensitivity","Double Vision", "Restlessness")
## Raw data plot for manuscript ##

# Prep data
raw_rpq_plot_df <- cbind.data.frame(factor(df$group, 
                                           levels = c("Military Control","Breacher", "Sniper")),
                                    mod_rpq)

# Reshape
raw_rpq_plot_long_df <- reshape(
  raw_rpq_plot_df,
  varying = list(names(raw_rpq_plot_df[2:17])),
  v.names = "response",
  timevar = "question",
  times = names(raw_rpq_plot_df[2:17]),  
  idvar = "id",          
  direction = "long"
)

# Change group column name
colnames(raw_rpq_plot_long_df)[1] <- "group"

# Create cumulative distribution function (CDF) data frame
raw_cdf_df <- raw_rpq_plot_long_df %>%
  group_by(question, group, response) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(question, group) %>%
  arrange(response) %>%
  mutate(
    total = sum(n),
    cum_prop = cumsum(n) / total * 100
  )

# Plot
raw_cdf_plot <- ggplot(raw_cdf_df, aes(x = response, y = cum_prop, color = group)) +
  geom_line() +
  geom_point() +
  facet_wrap(~question, scales = "free", ncol = 4) +
  scale_color_manual(values = c("orange", "blue", "black")) +
  scale_x_continuous(breaks = 0:5, name = "Response") +
  scale_y_continuous(name = "Cumulative Proportion (%)", limits = c(0, 100)) +
  theme_classic() +
  theme(legend.position = "top", legend.title = element_blank(),
        legend.text = element_text(face = "bold", size = 11),
        strip.background = element_blank(), strip.placement = "outside",
        strip.text = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"))  
raw_cdf_plot

# Save to file
ggsave("fig5.jpg", raw_cdf_plot, 
       width = 10, height = 10, units = "in",dpi = 600)

### Ordinal Logit Modelling ###

# List each RPQ3 score as an outcome in a datalist, dealing with NAs
dat_rpq <- list()
for (i in 1:16){
  dat_rpq[[i]] <- list(
    G = length(unique(df$group)),
    outcome = as.integer(mod_rpq[,i][complete.cases(mod_rpq[, i])] + 1),
    group = as.integer(factor(df$group[!(is.na(mod_rpq[, i]))], levels = c("Military Control",
                                                                           "Breacher", "Sniper"))),
    n = length(mod_rpq[,i][!(is.na(mod_rpq[, i]))]),
    K = 6)
}

dat_rpq[1:16]

# Compile model
ord_logit <- stan_model("ordinal_group_fixed.stan")

# Sample from data model
m_rpq <- list()
for (i in 1:16){
  m_rpq[[i]] <- sampling(ord_logit, data = dat_rpq[[i]], 
                          iter = 2000, 
                          chains = 4)}

# Check
m_rpq[[1]]

# Check model diagnostics
diagnostics <- list()
check <- list()
for (i in 1:length(m_rpq)){
diagnostics[[i]] <- util$extract_hmc_diagnostics(m_rpq[[i]])
check[[i]] <- util$check_all_hmc_diagnostics(diagnostics[[i]])}

## Extract samples across the list and evaluate group differences on cumulative
## probabilities. This only works for a 2 group situation.
m_rpq_posterior <- lapply(m_rpq, function(x) rstan::extract(x))

# Name the lists by the first 22 columns of the symp_comp_df3 dataframe
names(m_rpq_posterior) <- colnames(mod_rpq)

head(m_rpq_posterior[[14]])

m_rpq[[14]]

## Extract cumulative probabilities for each question, by group

# Function to derive the CDF
get_param_cdf_by_group <- function(gamma_draws, cutpoint_draws) {
  n_iter <- nrow(gamma_draws)
  G <- ncol(gamma_draws)
  K <- ncol(cutpoint_draws) + 1
  
  cdf_list <- vector("list", G)
  
  for (g in 1:G) {
    cdf_g <- matrix(NA, nrow = n_iter, ncol = K)
    
    for (s in 1:n_iter) {
      gamma_sg <- gamma_draws[s, g]
      cutpoints_s <- cutpoint_draws[s, ]
      
      # Apply logistic CDF for each cutpoint minus the group's gamma
      cdf_vals <- plogis(cutpoints_s - gamma_sg)
      
      cdf_g[s, 1:(K - 1)] <- cdf_vals
      cdf_g[s, K] <- 1  # The final category always has cumulative probability 1
    }
    
    cdf_list[[g]] <- cdf_g
  }
  
  return(cdf_list)
}

# Run function across list
cdf_list <- lapply(m_rpq_posterior, function(x) {
  get_param_cdf_by_group(x$gamma, x$cut_points)
})

cdf_list[[1]]

# Create group labels to label sublist elements
group_labels <- c("Military Control", "Breacher", "Sniper")

# Create data frames and add question and group info
for (i in seq_along(cdf_list)) {
  question_name <- names(cdf_list)[i]
  
  # Convert each group's matrix to data frame and label
  cdf_list[[i]] <- lapply(seq_along(cdf_list[[i]]), function(j) {
    df <- as.data.frame(cdf_list[[i]][[j]])
    df$question <- question_name
    df$group <- group_labels[j]
    return(df)
  })
}

# Take question level dataframes and combine into a single dataframe
cdf_df <- do.call(rbind, unlist(cdf_list, recursive = FALSE))
colnames(cdf_df)[c(1:6)] <- c(0, 1, 2, 3, 4, 5)

## Posterior Predictive Checks ##

# Pull the y_pred values out of the model posteriors
y_pred_list <- lapply(m_rpq_posterior, function(x) x$y_pred)

# Loop over questions
cdf_draws_list <- list()

for (q in seq_along(y_pred_list)) {
  question_name <- names(y_pred_list)[q]
  pred_mat <- y_pred_list[[q]]  # 4000 x n_obs for this question
  n_draws <- nrow(pred_mat)
  n_obs <- ncol(pred_mat)
  
  # Subset group vector to match participants used in this model
  group_subset <- df$group[!is.na(mod_rpq[[question_name]])]
  group_subset <- factor(group_subset, levels = c("Military Control", "Breacher", "Sniper"))
  
  # Defensive check
  stopifnot(length(group_subset) == n_obs)
  
  # Split column indices by group
  group_indices <- split(seq_along(group_subset), group_subset)
  
  # Store per-group results
  group_results <- list()
  
  for (g in names(group_indices)) {
    cols <- group_indices[[g]]
    submat <- pred_mat[, cols, drop = FALSE]
    
    # Compute cumulative proportions
    cdf_matrix <- apply(submat, 1, function(row) {
      tab <- table(factor(row, levels = 1:6))
      cumsum(as.numeric(tab)) / length(row) * 100
    })
    
    # Format to long
    cdf_df2 <- as.data.frame(t(cdf_matrix))
    colnames(cdf_df2) <- as.character(0:5)
    cdf_df2$draw <- 1:nrow(cdf_df2)
    cdf_df2$group <- g
    cdf_df2$question <- question_name
    
    cdf_long <- pivot_longer(
      cdf_df2,
      cols = all_of(as.character(0:5)),
      names_to = "response",
      values_to = "cum_prop"
    )
    
    group_results[[g]] <- cdf_long
  }
  
  # Combine group results
  cdf_draws_list[[q]] <- dplyr::bind_rows(group_results)
}

# Final combined PPC data
ppc_cdf_df <- do.call(rbind, cdf_draws_list)

# Subset to 100 draws for clarity
set.seed(123)
draw_subset <- sample(unique(ppc_cdf_df$draw), 100)
ppc_cdf_df_subset <- ppc_cdf_df[ppc_cdf_df$draw %in% draw_subset, ]

# Add overlay to existing plot
rpq_ppc_plot <- raw_cdf_plot +
  geom_line(data = ppc_cdf_df_subset,
            aes(x = as.numeric(response), y = cum_prop,
                group = interaction(draw, group), color = group),
            alpha = 0.05, linewidth = 0.3, inherit.aes = FALSE)
rpq_ppc_plot

# Second PPC without line draws
# Convert response to numeric just in case
ppc_cdf_df$response <- as.numeric(as.character(ppc_cdf_df$response))

# Summarize to mean and 90% interval
ppc_cdf_summary <- ppc_cdf_df %>%
  group_by(question, group, response) %>%
  summarise(
    lower = quantile(cum_prop, 0.05),
    upper = quantile(cum_prop, 0.95),
    mean  = mean(cum_prop),
    .groups = "drop"
  )

# Plot PPC2
group_colors <- c("Military Control" = "orange", 
                  "Breacher" = "blue", 
                  "Sniper" = "black")

rpq_ppc_summary_plot <- ggplot() +
  # Raw data dots
  geom_point(data = raw_cdf_df,
             aes(x = response, y = cum_prop, color = group),
             alpha = 0.7) +
  
  # 90% credible interval ribbon
  geom_ribbon(data = ppc_cdf_summary,
              aes(x = response, ymin = lower, ymax = upper, fill = group),
              alpha = 0.2, color = NA) +
  
  # Posterior mean line
  geom_line(data = ppc_cdf_summary,
            aes(x = response, y = mean, color = group),
            linewidth = 0.7) +
  
  facet_wrap(~question, scales = "free", ncol = 4) +
  
  # Harmonized scale
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  
  # One unified legend
  guides(fill = "none", color = guide_legend(override.aes = list(alpha = 1))) +
  
  scale_x_continuous(breaks = 0:5, name = "Response") +
  scale_y_continuous(name = "Cumulative Proportion (%)", limits = c(0, 100)) +
  theme_classic() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(face = "bold", size = 11),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  )
rpq_ppc_summary_plot 

## Table and results creation; comparing groups on their asymptomatic probability
## (probability of choosing zero) and also their probability of choosing three or higher,
## which is a surrogate of "high" symptoms

## Probability of choosing zero ##

# take cdf_df and create table-ready df
cdf_zero_df <- cdf_df[c(1, 7, 8)]
cdf_zero_df$`0` <- cdf_zero_df$`0` * 100
colnames(cdf_zero_df)[1] <- "value"

## Reshape wide ##

# Add draw ID
cdf_zero_df  <- cdf_zero_df  %>%
  group_by(group, question) %>%
  mutate(draw = row_number()) %>%
  ungroup()

# Pivot to wide format
cdf_zero_df_wide <- cdf_zero_df  %>%
  select(group, draw, question, value) %>%
  pivot_wider(names_from = question, values_from = value)

# Re order group
cdf_zero_df_wide$group <- factor(cdf_zero_df_wide$group, 
                              levels = c("Military Control", "Breacher", "Sniper"))

# Make table
rpq_zero_tbl <-
  tbl_summary(cdf_zero_df_wide[-c(2)],
              by = group,
              missing = "no",
              digits = all_continuous() ~ 1,
              statistic = all_continuous()~c("{mean} ({p5} - {p95})")) %>%
  modify_header(label = "**Symptom**", all_stat_cols() ~ "**{level}**") %>%
  modify_table_styling(columns = everything(), footnote = NA) %>%
  bold_labels() 

# Creating Contrast df's 
bvm_contrast <- cdf_zero_df_wide[-c(1,2)][cdf_zero_df_wide$group=="Military Control",] -
  cdf_zero_df_wide[-c(1,2)][cdf_zero_df_wide$group== "Breacher",]
svm_contrast <- cdf_zero_df_wide[-c(1,2)][cdf_zero_df_wide$group=="Military Control",] -
  cdf_zero_df_wide[-c(1,2)][cdf_zero_df_wide$group=="Sniper",]
bvs_contrast <- cdf_zero_df_wide[-c(1,2)][cdf_zero_df_wide$group=="Breacher",] -
  cdf_zero_df_wide[-c(1,2)][cdf_zero_df_wide$group=="Sniper",]

# Create contrast tables

# Function to get posterior prob > 0
prob_func <- function(x) {
  out <- mean(x > 0) * 100
  return(out)
}

bvm_tbl <- 
  tbl_summary(bvm_contrast,
              missing = "no",
              digits = all_continuous() ~ 1,
              statistic = all_continuous()~c("{mean} ({p5} - {p95}; pprob = {pprob_func})")) %>%
  bold_labels() %>%
  modify_header(label = "**Symptom**", stat_0 = "**Military Control – Breacher**") %>%
  modify_table_styling(columns = everything(), footnote = NA)

svm_tbl <- 
  tbl_summary(svm_contrast,
              missing = "no",
              digits = all_continuous() ~ 1,
              statistic = all_continuous()~c("{mean} ({p5} - {p95}; pprob = {pprob_func})")) %>%
  bold_labels() %>%
  modify_header(label = "**Symptom**", stat_0 = "**Military Control – Sniper**") %>%
  modify_table_styling(columns = everything(), footnote = NA)

bvs_tbl <- 
  tbl_summary(bvs_contrast,
              missing = "no",
              digits = all_continuous() ~ 1,
              statistic = all_continuous()~c("{mean} ({p5} - {p95}; pprob = {pprob_func})")) %>%
  bold_labels() %>%
  modify_header(label = "**Symptom**", stat_0 = "**Breacher – Sniper**") %>%
  modify_table_styling(columns = everything(), footnote = NA)

contrast_tbl <- 
  tbl_merge(list(bvm_tbl, svm_tbl, bvs_tbl))

#merging raw values and posterior estimates
rpq0_tbl <-
  tbl_merge(list(rpq_zero_tbl, contrast_tbl),
            tab_spanner = c("**Estimated Probability of Choosing Zero**", 
                            "**Differences (%)**")) %>%
  as_gt()%>%
  tab_source_note(source_note = md('The table displays the estimated percentage of 
  participants who reported no symptoms for each group, along
  with the posterior estimates of the differences in percentages between groups.'))%>%
  tab_source_note(source_note = md('Values displayed as the estimated mean (90% Credible Interval), 
                                   with differences presented as the estimated mean (90% CI; pprob = 
                                   percentage of the posterior density above 0%).'))%>%
  tab_source_note(source_note = md('Summary statistics are derived from 2000 posterior samples.'))%>%
  
  gt::gtsave(filename = 'sup_tbl2.html')


## Probability of choosing 3 or higher ##

# take cdf_df and create table-ready df
# in this instance, we take the cut off value of 2, and subtract it from 1
# Therefore =  1 - 2nd cut_point
cdf_3plus_df <- cdf_df[c(3, 7, 8)]
cdf_3plus_df$`2` <- 1 - cdf_3plus_df$`2` # 1 - 2nd cut point
cdf_3plus_df$`2` <- cdf_3plus_df$`2` * 100 # values as percentages
colnames(cdf_3plus_df)[1] <- "value"

# Reshape wide

# Add draw ID
cdf_3plus_df  <- cdf_3plus_df  %>%
  group_by(group, question) %>%
  mutate(draw = row_number()) %>%
  ungroup()

# Pivot to wide format
cdf_3plus_df_wide <- cdf_3plus_df  %>%
  select(group, draw, question, value) %>%
  pivot_wider(names_from = question, values_from = value)

# Re order group
cdf_3plus_df_wide$group <- factor(cdf_3plus_df_wide$group, 
                                 levels = c("Military Control", "Breacher", "Sniper"))

# Make table
rpq_3plus_tbl <-
  tbl_summary(cdf_3plus_df_wide[-c(2)],
              by = group,
              missing = "no",
              digits = all_continuous() ~ 1,
              statistic = all_continuous()~c("{mean} ({p5} - {p95})")) %>%
  modify_header(label = "**Symptom**", all_stat_cols() ~ "**{level}**") %>%
  modify_table_styling(columns = everything(), footnote = NA) %>%
  bold_labels() 

# Creating Contrast df's 
bvm3_contrast <- cdf_3plus_df_wide[-c(1,2)][cdf_3plus_df_wide$group == "Breacher",] -
  cdf_3plus_df_wide[-c(1,2)][cdf_3plus_df_wide$group == "Military Control",]
svm3_contrast <- cdf_3plus_df_wide[-c(1,2)][cdf_3plus_df_wide$group=="Sniper",] -
  cdf_3plus_df_wide[-c(1,2)][cdf_3plus_df_wide$group == "Military Control",]
bvs3_contrast <- cdf_3plus_df_wide[-c(1,2)][cdf_3plus_df_wide$group == "Breacher",] -
  cdf_3plus_df_wide[-c(1,2)][cdf_3plus_df_wide$group == "Sniper",]

# Create contrast tables

# Function to get posterior prob > 0
prob_func <- function(x) {
  out <- mean(x > 0) * 100
  return(out)
}

bvm3_tbl <- 
  tbl_summary(bvm3_contrast,
              missing = "no",
              digits = all_continuous() ~ 1,
              statistic = all_continuous()~c("{mean} ({p5} - {p95}; pprob = {pprob_func})")) %>%
  bold_labels() %>%
  modify_header(label = "**Symptom**", stat_0 = "**Breacher - Military Control**") %>%
  modify_table_styling(columns = everything(), footnote = NA)

svm3_tbl <- 
  tbl_summary(svm3_contrast,
              missing = "no",
              digits = all_continuous() ~ 1,
              statistic = all_continuous()~c("{mean} ({p5} - {p95}; pprob = {pprob_func})")) %>%
  bold_labels() %>%
  modify_header(label = "**Symptom**", stat_0 = "**Sniper - Military Control**") %>%
  modify_table_styling(columns = everything(), footnote = NA)

bvs3_tbl <- 
  tbl_summary(bvs3_contrast,
              missing = "no",
              digits = all_continuous() ~ 1,
              statistic = all_continuous()~c("{mean} ({p5} - {p95}; pprob = {pprob_func})")) %>%
  bold_labels() %>%
  modify_header(label = "**Symptom**", stat_0 = "**Breacher – Sniper**") %>%
  modify_table_styling(columns = everything(), footnote = NA)

contrast3_tbl <- 
  tbl_merge(list(bvm3_tbl, svm3_tbl, bvs3_tbl))

#merging raw values and posterior estimates
rpq3_tbl <-
  tbl_merge(list(rpq_3plus_tbl, contrast3_tbl),
            tab_spanner = c("**Estimated Probability of Choosing 3 or Higher**", 
                            "**Differences (%)**")) %>%
  as_gt()%>%
  tab_source_note(source_note = md('The table displays the estimated percentage of 
  participants who reported a symptom score of 3 or higher for each group, along
  with the posterior estimates of the differences in percentages between groups.'))%>%
  tab_source_note(source_note = md('Values displayed as the estimated mean (90% Credible Interval), 
                                   with differences presented as the estimated mean (90% CI; pprob = 
                                   percentage of the posterior density above 0%).'))%>%
  tab_source_note(source_note = md('Summary statistics are derived from 2000 posterior samples.'))%>%
  
  gt::gtsave(filename = 'sup_tbl3.html')

























