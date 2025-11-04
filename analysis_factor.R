## This script is used to model the latent cognitive variable
## that is used downstream to compare breachers/snipers with their military
## control counterparts

#Libraries Used
library(tidyverse)
library(magrittr)
library(modelr)
library(ggplot2)
library(cowplot)
library(rstan)
library(rstanarm)
library(rethinking)  
library(tidybayes)
library(ggrepel)
library(RColorBrewer)
library(gganimate)
library(loo)
library(StanHeaders)

# Functions
z_score <- function(x) {
  out <- (x - mean(x,na.rm=TRUE))/sd(x,na.rm = TRUE)}
mod_z <- function(x) {
  out <- (x - median(x,na.rm=TRUE))/mad(x,na.rm = TRUE)}

# For speed with mcmc
options(mc.cores = parallel::detectCores())

# Load file (Not available for public release)
df <- read.csv('data_df.csv', stringsAsFactors = FALSE)
colnames(df)

# Isolate just the cognitive variables for modelling
cog_df <- df[c(17:22)]

# Create the elements of the data list
n <- nrow(cog_df)
m <- ncol(cog_df)
n_miss_cog <- length(cog_df[is.na(cog_df)])

# Isolate missing from observed
Y_cog_df <-data.frame(apply(cog_df, 2, z_score))
na_values_cog <- is.na(Y_cog_df)
y_obs_cog <- Y_cog_df[!na_values_cog]
n_obs_cog <- length(y_obs_cog)
row_miss_cog <- rep(1:n, times = m)[na_values_cog]
col_miss_cog <- rep(1:m, each = n)[na_values_cog]
miss_check_cog <- cbind(row_miss_cog, col_miss_cog)

# Identifies index where observed values are
row_obs_cog <- rep(1:n, times = m)[!na_values_cog]
col_obs_cog <- rep(1:m, each = n)[!na_values_cog]

dat_cog <- list(row_obs = row_obs_cog,
                k = 1,
                m = m,
                n = n,
                col_obs = col_obs_cog,
                row_miss = row_miss_cog,
                col_miss = col_miss_cog,
                n_miss = n_miss_cog,
                n_obs = n_obs_cog,
                y_obs = y_obs_cog,
                y = as.matrix(Y_cog_df))

str(dat_cog)

# Compile models
factor_impute_mcar <- stan_model("factor_analysis_impute_mcar.stan")
mod_prior <- stan_model("factor_analysis_prior.stan")

m_cog <- sampling(object = factor_impute_mcar, data = dat_cog, 
                  chains = 4, cores = 10,
                  iter = 3000, init_r = .2) # add some init


# Sample the prior to evaluate prior choices
# Create fake matrix 
y_prior <- matrix(0, nrow = n, ncol = m)

dat_prior <- list(
  k = 1,
  m = m,
  n = n,
  n_obs = n * m,
  y = y_prior)

str(dat_prior)

# Sample from prior
m_cog_prior <- sampling(object = mod_prior, data = dat_prior, chains = 4, 
                        cores = 10,
                        iter = 3000, 
                        init_r = .2) # add some init
m_cog_prior

# Extract the posterior
post_cog <- data.frame(extract.samples(m_cog))
colnames(post_cog)
precis(post_cog, 2)

# Compare post factor loadings
# Create a df for factor loadings
post_L_cog <- post_cog[grepl("L", colnames(post_cog)) &
                         !(grepl("Sigma", colnames(post_cog)))]
# Sigma leftover, just remove.
post_L_cog$sigma_L <- NULL

# Create separate dataframes, and then combine
f1_cog <- post_L_cog[c(1:m)]
colnames(cog_df)
colnames(f1_cog) <- c("4-choice RT task (Mean RT Correct)",
                      "Delayed Matching-to-Sample Task (% Accurate)",
                      "1-Back", "2-Back", "3-Back", "Stroop (RT Difference)")

# Transform to long
f1_cog_long <- gather(f1_cog, key = "Cognitive Measure", value = "Loading score (SD Units)")

f1_cog_long$`Cognitive Measure` <- 
  factor(f1_cog_long$`Cognitive Measure`,
         levels = c("Stroop (RT Difference)",
                    "Delayed Matching-to-Sample Task (% Accurate)",
                    "4-choice RT task (Mean RT Correct)",
                    "1-Back","2-Back","3-Back"))

# Add a factor column
f1_cog_long$factor <- "Factor 1"

# Prior
# Extract the prior
prior_cog <- data.frame(extract.samples(m_cog_prior))
precis(prior_cog, 2)

# Compare prior factor loadings
# Create a df for factor loadings
prior_L_cog <- prior_cog[grepl("L", colnames(prior_cog)) &
                           !(grepl("Sigma", colnames(prior_cog)))]
# Sigma leftover, just remove.
prior_L_cog$sigma_L <- NULL

# Create separate dataframes, and then combine
f1_cog_prior <- prior_L_cog[c(1:m)]
colnames(f1_cog_prior) <- c("4-choice RT task (Mean RT Correct)",
                            "Delayed Matching-to-Sample Task (% Accurate)",
                            "1-Back", "2-Back", "3-Back", "Stroop (RT Difference)")

# Transform to long
f1_cog_long_prior <- gather(f1_cog_prior, key = "Cognitive Measure", value = "Loading score (SD Units)")
# Add factor #
f1_cog_long_prior$factor <- "Factor 1"

f1_cog_long_prior$`Cognitive Measure` <-
  factor(f1_cog_long_prior$`Cognitive Measure`,
         levels = c("Stroop (RT Difference)", "Delayed Matching-to-Sample Task (% Accurate)",
                    "4-choice RT task (Mean RT Correct)", "1-Back", "2-Back", "3-Back"))

# Combine
f1_cog_long$samples <- "Posterior"
f1_cog_long_prior$samples <- "Prior"

plot_cog <- rbind.data.frame(f1_cog_long, f1_cog_long_prior)

# Combined plot of prior and posterior samples
theme_set(theme_tidybayes() + panel_border())
plot_cog_loadings <- plot_cog %>%
  ggplot(aes(y = fct_rev(`Cognitive Measure`),
             x = `Loading score (SD Units)`,
             fill = samples)) +
  geom_vline(xintercept = 0, alpha = 0.8, linetype = 2) +
  theme(strip.text.x = element_text(face = "bold", size = 10),
        axis.title.x = element_text(size = 13, face = "bold"),
        axis.title.y = element_text(size = 13, face = "bold"),
        axis.text.y = element_text(face = "bold", size = 10),
        legend.position = "None",
        axis.text.x = element_text(face = "bold", size = 10)) +
  ylab("") + xlim (-2,2) +
  scale_fill_manual(values = c("pink","grey90")) +
  stat_halfeye(.width = c(0.70, 0.90), alpha = 0.7)
plot_cog_loadings
theme_set(theme_tidybayes() + panel_border())

ggsave("fig2.jpg", plot_cog_loadings, dpi = 600, width = 11)

# Save latent model files for downstream student_t modelling
write.csv(post_cog, "post_cog.csv", row.names = FALSE)


