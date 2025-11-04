### Code to run the Student T models based on the latent weights estimated in
### analysis_factor.R, by group

#libraries
library(rstan)
library(rstanarm)
library(rethinking)  
library(StanHeaders)

# Set options for rstan
# For speed with mcmc
options(mc.cores = parallel::detectCores())

## Load Files

# Load latent models derived from analysis_factor.R
post_cog <- read.csv("post_cog.csv", stringsAsFactors = FALSE)

# Load original data (Not available for public release)
df <- read.csv('data_df.csv', stringsAsFactors = FALSE)
table(df$group)

# Downstream Modelling of dag1 in dags.R using the individual latent factor scores
# to represent cognition and mental health

## Evaluate and create prior models to underlay on the posterior plots##

# Create data for prior (total effects)
dat_tot_prior <- list(group = as.integer(as.factor(df$group)),
                      G = 3,
                      outcome = rep(0, nrow(df)),
                      n = nrow(df))
str(dat_tot_prior)

#create data for prior (direct effects)
dat_dir_prior <- list(group = as.integer(as.factor(df$group)),
                      G = 3,
                      outcome = rep(0, nrow(df)),
                      ord_pred = rep(0, nrow(df)),
                      K = 3, # Levels in the BSI/PTSD variables
                      age = rep(0, nrow(df)),
                      n = nrow(df))
str(dat_dir_prior)

# Compile models (including zero inflated)
st_tot_prior <- stan_model("student_total_prior.stan")
st_dir_prior <- stan_model("student_direct_prior.stan")

# Sample prior models
m_tot_cog_prior <- sampling(st_tot_prior, data = dat_tot_prior, 
                            chains = 4, iter = 1000)
m_dir_cog_prior <- sampling(st_dir_prior, data = dat_dir_prior, 
                            chains = 4, iter = 1000)

## Modelling the data ##

# Derive the z parameter estimates (the weighted latent scores) from the cognitive and mh models
post_z_cog <- post_cog[grepl("z", colnames(post_cog))]

## Full cognitive total effects model with 100 draws from the latent model posterior##
# Isolate z scores
n_subsets <- 100
post_size <- nrow(post_z_cog)
subset_rows <- sample(post_size, n_subsets, replace = F )
z1_subset_tot_cog <- post_z_cog[subset_rows, ]

z1_samples_tot_cog <- list()
for (i in 1:100) {
  z1_samples_tot_cog[[i]] <- as.vector(t(z1_subset_tot_cog[i, ]))
}

#check
z1_samples_tot_cog[[5]]


#set up stan data list
stan_data_tot_cog <- list()
for (i in 1:100){
  n <- ncol(post_z_cog) #subjects
  G <- 3
  group <- as.integer(as.factor(df$group))#groups
  outcome <- z_score(z1_samples_tot_cog[[i]])
  stan_data_tot_cog[[i]] <- list(n = n, G = G, group = group, outcome = outcome)
}

#check
stan_data_tot_cog[[5]]

# Compile model
st_tot <- stan_model("student_total.stan")

#Run the loop to generate 100 samples and run each of these 100 models for 500 iterations
draws_tot_cog <- list()
for (i in 1:100) {
  draws_tot_cog[[i]] <- sampling(st_tot, 
                                 data = stan_data_tot_cog[[i]],  
                                 chains = 4, 
                                 iter = 500)
}

post_draws_tot_cog <- lapply(draws_tot_cog, as.data.frame(rstan::extract))
post_draws_tot_cog_df <- do.call(rbind, post_draws_tot_cog)
precis(post_draws_tot_cog_df, 2)

## Check model performance ##

# R hat values 
rhat_df <- list()
for (i in 1:100){
  df_prep <- as.data.frame(summary(draws_tot_cog[[i]])$summary)
  rhat_df[[i]] <- df_prep[c(9, 10)]
}

# Put together
rhat_tot_cog_df <- do.call(cbind, rhat_df)

# Make random draws
draws <- sample(1:nrow(post_draws_tot_cog_df), size = 1000)
colnames(post_draws_tot_cog_df)

# Pairs plot of parameters
post_pairs_tot_cog <- post_draws_tot_cog_df[draws, c(2:4)]
pairs(post_pairs_tot_cog)

### Direct effects model with adjustements, with 100 draws from ###
## the latent model posterior ###

# Isolate z scores
n_subsets <- 100
post_size <- nrow(post_cog)
subset_rows <- sample(post_size, n_subsets, replace = F )
z1_subset_dir_cog <- post_z_cog[subset_rows, ]

z1_samples_dir_cog <- list()
for (i in 1:100) {
  z1_samples_dir_cog[[i]] <- as.vector(t(z1_subset_dir_cog[i, ]))
}

# Check
z1_samples_dir_cog[[1]]

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

# Set up stan data lists
stan_data_dir_cog_ptsd <- list()
for (i in 1:100){
  n <- ncol(post_z_cog) #subjects
  G <- 3
  group <- as.integer(as.factor(df$group)) #groups
  ord_pred <- as.integer(df$PTSD_cat) # ordinal PTSD scores
  age <- z_score(df$Age)
  outcome <- z_score(z1_samples_dir_cog[[i]])
  stan_data_dir_cog_ptsd[[i]] <- list(n = n, G = G, group = group, K = 3,
                                  ord_pred = ord_pred,
                                 age = age, outcome = outcome)
}


stan_data_dir_cog_bsis <- list()
for (i in 1:100){
  n <- ncol(post_z_cog) #subjects
  G <- 3
  group <- as.integer(as.factor(df$group)) #groups
  ord_pred <- as.integer(df$bsis_cat) # ordinal PTSD scores
  age <- z_score(df$Age)
  outcome <- z_score(z1_samples_dir_cog[[i]])
  stan_data_dir_cog_bsis[[i]] <- list(n = n, G = G, group = group, K = 3,
                                      ord_pred = ord_pred,
                                      age = age, outcome = outcome)
}

stan_data_dir_cog_bsia <- list()
for (i in 1:100){
  n <- ncol(post_z_cog) #subjects
  G <- 3
  group <- as.integer(as.factor(df$group)) #groups
  ord_pred <- as.integer(df$bsia_cat) # ordinal PTSD scores
  age <- z_score(df$Age)
  outcome <- z_score(z1_samples_dir_cog[[i]])
  stan_data_dir_cog_bsia[[i]] <- list(n = n, G = G, group = group, K = 3,
                                      ord_pred = ord_pred,
                                      age = age, outcome = outcome)
}

stan_data_dir_cog_bsid <- list()
for (i in 1:100){
  n <- ncol(post_z_cog) #subjects
  G <- 3
  group <- as.integer(as.factor(df$group)) #groups
  ord_pred <- as.integer(df$bsid_cat) # ordinal PTSD scores
  age <- z_score(df$Age)
  outcome <- z_score(z1_samples_dir_cog[[i]])
  stan_data_dir_cog_bsid[[i]] <- list(n = n, G = G, group = group, K = 3,
                                      ord_pred = ord_pred,
                                      age = age, outcome = outcome)
}

# Compile
st_dir <- stan_model("student_direct.stan")

# Sample from models
#Run the loop to generate 200 samples and run each of these 200 models for 1000 iterations
draws_dir_cog_ptsd <- list()
for (i in 1:100) {
  draws_dir_cog_ptsd[[i]] <- sampling(st_dir, 
                                 data = stan_data_dir_cog_ptsd[[i]],
                                 chains = 4, 
                                 iter = 500)
}

draws_dir_cog_bsis <- list()
for (i in 1:100) {
  draws_dir_cog_bsis[[i]] <- sampling(st_dir, 
                                 data = stan_data_dir_cog_bsis[[i]],
                                 chains = 4, 
                                 iter = 500)
}

draws_dir_cog_bsia <- list()
for (i in 1:100) {
  draws_dir_cog_bsia[[i]] <- sampling(st_dir, 
                                      data = stan_data_dir_cog_bsia[[i]],
                                      chains = 4, 
                                      iter = 500)
}

draws_dir_cog_bsid <- list()
for (i in 1:100) {
  draws_dir_cog_bsid[[i]] <- sampling(st_dir, 
                                      data = stan_data_dir_cog_bsid[[i]],
                                      chains = 4, 
                                      iter = 500)
}

### Posterior summaries and checks ###

## Extract posteriors ##

# PTSD #
post_draws_dir_cog_ptsd <- lapply(draws_dir_cog_ptsd, as.data.frame(extract.samples))
post_draws_dir_cog_ptsd_df <- do.call(rbind, post_draws_dir_cog_ptsd)
precis(post_draws_dir_cog_ptsd_df, 2)

# BSI-S #
post_draws_dir_cog_bsis <- lapply(draws_dir_cog_bsis, as.data.frame(extract.samples))
post_draws_dir_cog_bsis_df <- do.call(rbind, post_draws_dir_cog_bsis)
precis(post_draws_dir_cog_bsis_df, 2)

# BSI-A #
post_draws_dir_cog_bsia <- lapply(draws_dir_cog_bsia, as.data.frame(extract.samples))
post_draws_dir_cog_bsia_df <- do.call(rbind, post_draws_dir_cog_bsia)
precis(post_draws_dir_cog_bsia_df, 2)

# BSI-D #
post_draws_dir_cog_bsid <- lapply(draws_dir_cog_bsid, as.data.frame(extract.samples))
post_draws_dir_cog_bsid_df <- do.call(rbind, post_draws_dir_cog_bsid)
precis(post_draws_dir_cog_bsid_df, 2)


## Check model performance ##

# PTSD #

# r hat values 
rhat_df <- list()
for (i in 1:100){
  df_prep <- as.data.frame(summary(draws_dir_cog_ptsd[[i]])$summary)
  rhat_df[[i]] <- df_prep[c(9, 10)]
}

# Put together
rhat_dir_cog_ptsd_df <- do.call(cbind, rhat_df)

# Make random draws
draws <- sample(1:nrow(post_draws_dir_cog_ptsd_df), size = 1000)
colnames(post_draws_dir_cog_ptsd_df)

# Pairs plot of parameters
post_pairs_dir_cog_ptsd <- post_draws_dir_cog_ptsd_df[draws, c(2:4)] # Change to see others
pairs(post_pairs_dir_cog_ptsd)

# BSI-S #

# r hat values 
rhat_df <- list()
for (i in 1:100){
  df_prep <- as.data.frame(summary(draws_dir_cog_bsis[[i]])$summary)
  rhat_df[[i]] <- df_prep[c(9, 10)]
}

# Put together
rhat_dir_cog_bsis_df <- do.call(cbind, rhat_df)

# Make random draws
draws <- sample(1:nrow(post_draws_dir_cog_bsis_df), size = 1000)
colnames(post_draws_dir_cog_bsis_df)

# Pairs plot of parameters
post_pairs_dir_cog_bsis <- post_draws_dir_cog_bsis_df[draws, c(2:4)] # Change to see others
pairs(post_pairs_dir_cog_bsis)

# BSI-A #

# r hat values 
rhat_df <- list()
for (i in 1:100){
  df_prep <- as.data.frame(summary(draws_dir_cog_bsia[[i]])$summary)
  rhat_df[[i]] <- df_prep[c(9, 10)]
}

# Put together
rhat_dir_cog_bsia_df <- do.call(cbind, rhat_df)

# Make random draws
draws <- sample(1:nrow(post_draws_dir_cog_bsia_df), size = 1000)
colnames(post_draws_dir_cog_bsia_df)

# Pairs plot of parameters
post_pairs_dir_cog_bsia <- post_draws_dir_cog_bsia_df[draws, c(2:4)] # Change to see others
pairs(post_pairs_dir_cog_bsia)

# BSI-D #

# r hat values 
rhat_df <- list()
for (i in 1:100){
  df_prep <- as.data.frame(summary(draws_dir_cog_bsid[[i]])$summary)
  rhat_df[[i]] <- df_prep[c(9, 10)]
}

# Put together
rhat_dir_cog_bsid_df <- do.call(cbind, rhat_df)

# Make random draws
draws <- sample(1:nrow(post_draws_dir_cog_bsid_df), size = 1000)
colnames(post_draws_dir_cog_bsid_df)

# Pairs plot of parameters
post_pairs_dir_cog_bsid <- post_draws_dir_cog_bsid_df[draws, c(2:4)] # Change to see others
pairs(post_pairs_dir_cog_bsid)


#### Plotting and Results reporting ####

## Generate a sample of 1000 rows from the posterior across ##
## all columns ##

# Total cog model
colnames(post_draws_tot_cog_df)
sample_tot_cog <- sample(1:nrow(post_draws_tot_cog_df),1000)
post_cog_tot_df <- post_draws_tot_cog_df[sample_tot_cog, c(2:4)]

# Direct cog model ptsd
colnames(post_draws_dir_cog_ptsd_df)
sample_dir_cog_ptsd <- sample(1:nrow(post_draws_dir_cog_ptsd_df),1000)
post_cog_dir_ptsd_df <- post_draws_dir_cog_ptsd_df[sample_dir_cog_ptsd, c(2:4)]

# Direct cog model bsis
colnames(post_draws_dir_cog_bsis_df)
sample_dir_cog_bsis <- sample(1:nrow(post_draws_dir_cog_bsis_df),1000)
post_cog_dir_bsis_df <- post_draws_dir_cog_bsis_df[sample_dir_cog_bsis, c(2:4)]

# Direct cog model bsia
colnames(post_draws_dir_cog_bsia_df)
sample_dir_cog_bsia <- sample(1:nrow(post_draws_dir_cog_bsia_df),1000)
post_cog_dir_bsia_df <- post_draws_dir_cog_bsia_df[sample_dir_cog_bsia, c(2:4)]

# Direct cog model bsid
colnames(post_draws_dir_cog_bsid_df)
sample_dir_cog_bsid <- sample(1:nrow(post_draws_dir_cog_bsid_df),1000)
post_cog_dir_bsid_df <- post_draws_dir_cog_bsid_df[sample_dir_cog_bsid, c(2:4)]

# Name columns for each
colnames(post_cog_tot_df) <- c("Breachers", "Military Controls", "Snipers")
colnames(post_cog_dir_ptsd_df) <- c("Breachers", "Military Controls", "Snipers")
colnames(post_cog_dir_bsis_df) <- c("Breachers", "Military Controls", "Snipers")
colnames(post_cog_dir_bsia_df) <- c("Breachers", "Military Controls", "Snipers")
colnames(post_cog_dir_bsid_df) <- c("Breachers", "Military Controls", "Snipers")

# Add contrasts #
# Total Cog Model
post_cog_tot_df$contrast_bc <- post_cog_tot_df$Breachers -
  post_cog_tot_df$`Military Controls`
post_cog_tot_df$contrast_bs <- post_cog_tot_df$Breachers -
  post_cog_tot_df$Snipers
post_cog_tot_df$contrast_sc <- post_cog_tot_df$Snipers -
  post_cog_tot_df$`Military Controls`

# Direct Cog Model ptsd
post_cog_dir_ptsd_df$contrast_bc <- post_cog_dir_ptsd_df$Breachers -
  post_cog_dir_ptsd_df$`Military Controls`
post_cog_dir_ptsd_df$contrast_bs <- post_cog_dir_ptsd_df$Breachers -
  post_cog_dir_ptsd_df$Snipers
post_cog_dir_ptsd_df$contrast_sc <- post_cog_dir_ptsd_df$Snipers -
  post_cog_dir_ptsd_df$`Military Controls`

# Direct Cog Model bsis
post_cog_dir_bsis_df$contrast_bc <- post_cog_dir_bsis_df$Breachers -
  post_cog_dir_bsis_df$`Military Controls`
post_cog_dir_bsis_df$contrast_bs <- post_cog_dir_bsis_df$Breachers -
  post_cog_dir_bsis_df$Snipers
post_cog_dir_bsis_df$contrast_sc <- post_cog_dir_bsis_df$Snipers -
  post_cog_dir_bsis_df$`Military Controls`

# Direct Cog Model bsia
post_cog_dir_bsia_df$contrast_bc <- post_cog_dir_bsia_df$Breachers -
  post_cog_dir_bsia_df$`Military Controls`
post_cog_dir_bsia_df$contrast_bs <- post_cog_dir_bsia_df$Breachers -
  post_cog_dir_bsia_df$Snipers
post_cog_dir_bsia_df$contrast_sc <- post_cog_dir_bsia_df$Snipers -
  post_cog_dir_bsia_df$`Military Controls`

# Direct Cog Model bsid
post_cog_dir_bsid_df$contrast_bc <- post_cog_dir_bsid_df$Breachers -
  post_cog_dir_bsid_df$`Military Controls`
post_cog_dir_bsid_df$contrast_bs <- post_cog_dir_bsid_df$Breachers -
  post_cog_dir_bsid_df$Snipers
post_cog_dir_bsid_df$contrast_sc <- post_cog_dir_bsid_df$Snipers -
  post_cog_dir_bsid_df$`Military Controls`

## Create prior draws to overlay on graphs

# Extract priors
post_cog_tot_prior <- data.frame(rstan::extract(m_tot_cog_prior))
post_cog_dir_prior <- data.frame(rstan::extract(m_dir_cog_prior))

# Match posterior columns
post_cog_tot_prior <- post_cog_tot_prior[c(2:4)]
post_cog_dir_prior <- post_cog_dir_prior[c(2:4)]

# Create contrasts for prior
post_cog_tot_prior$contrast_bc <- post_cog_tot_prior$g.1 -
  post_cog_tot_prior$g.2
post_cog_tot_prior$contrast_bs <- post_cog_tot_prior$g.1 -
  post_cog_tot_prior$g.3
post_cog_tot_prior$contrast_sc <- post_cog_tot_prior$g.3 -
  post_cog_tot_prior$g.2
post_cog_dir_prior$contrast_bc <- post_cog_dir_prior$g.1 -
  post_cog_dir_prior$g.2
post_cog_dir_prior$contrast_bs <- post_cog_dir_prior$g.1 -
  post_cog_dir_prior$g.3
post_cog_dir_prior$contrast_sc <- post_cog_dir_prior$g.3 -
  post_cog_dir_prior$g.2

# Name columns
colnames(post_cog_tot_prior) <- colnames(post_cog_tot_df)
colnames(post_cog_dir_prior) <- colnames(post_cog_tot_df)

###############PLOTS################

## Cog Total Plot ##
pdf('fig3.pdf', width = 9, height = 4.5, pointsize = 10)
par(mfrow = c(1, 2))

## A. Unadjusted Group Estimates
dens(post_cog_tot_df$Breachers, xlim = c(-1.5, 1.5),
     lwd = 3, col = "darkorange", 
     xlab = "Latent Neurocognitive Scores (Unadjusted, SD Units)", 
     ylim = c(0, 4), font.lab = 2, 
     cex.lab = 0.9,  frame = FALSE)
dens(post_cog_tot_df$Snipers, lwd = 3, col = "blue3", add = TRUE)
dens(post_cog_tot_df$`Military Controls`, lwd = 3, col = "forestgreen", add = TRUE)
dens(post_cog_tot_prior$Breachers, lwd = 1, lty = 6, col = 'gray80', add = TRUE)
dens(post_cog_tot_prior$Snipers, lwd = 1, lty = 6, col = 'gray60', add = TRUE)
dens(post_cog_tot_prior$`Military Controls`, lwd = 1, lty =
       6, col = 'gray40', add = TRUE)
legend("topright", box.lty = 0, cex = 0.9,
       legend = c("Breachers", "Snipers", "Military Controls"),
       fill = c("darkorange", "blue3", "forestgreen"))
text(-1.17, 0.34, "Breachers Prior", cex = 0.7, col = 'gray80')
text(-1.17, 0.5, "Snipers Prior", cex = 0.7, col = 'gray60')
text(-1.17, 0.17, "Military Controls Prior", cex = 0.7, col = 'gray40')
text(-1.5, 4, "A", cex = 0.9, font = 2)

## B. Unadjusted Contrasts (all on same plot)
plot(density(as.numeric(post_cog_tot_df$contrast_bc), adjust = 0.5),
     xlim = c(-1.5, 1.5), ylim = c(0, 4),
     lwd = 3, col = "firebrick", 
     main = "",
     font.lab = 2, cex.lab = 0.9, 
     xlab = "Difference in Latent Neurocognitive Scores (Unadjusted, SD Units)", 
     ylab = "", frame = FALSE)
lines(density(as.numeric(post_cog_tot_df$contrast_sc), adjust = 0.5),
      lwd = 3, col = "steelblue4")
lines(density(as.numeric(post_cog_tot_df$contrast_bs), adjust = 0.5),
      lwd = 3, col = "goldenrod3")
abline(v = 0, lty = 3, col = "gray40")
legend("topright", box.lty = 0, cex = 0.9,
       legend = c("Breachers - Military Controls", "Snipers - Military Controls", "Breachers - Snipers"),
       col = c("firebrick", "steelblue4", "goldenrod3"), lwd = 3)
text(-1.5, 4, "B", cex = 0.9, font = 2)

fig3 <- recordPlot()
dev.off()

## Cog Total Plot total with all adjustments##
pdf('sup_fig1.pdf', width = 9, height = 13, pointsize = 12)
par(mfrow = c(5, 2))

## A. Unadjusted Group Estimates
dens(post_cog_tot_df$Breachers, xlim = c(-1.5, 1.5),
     lwd = 3, col = "darkorange", 
     xlab = "Latent Neurocognitive Scores (Unadjusted, SD Units)", 
     ylim = c(0, 3), font.lab = 2, 
     cex.lab = 0.9,  frame = FALSE)
dens(post_cog_tot_df$Snipers, lwd = 3, col = "blue3", add = TRUE)
dens(post_cog_tot_df$`Military Controls`, lwd = 3, col = "forestgreen", add = TRUE)
dens(post_cog_tot_prior$Breachers, lwd = 1, lty = 6, col = 'gray80', add = TRUE)
dens(post_cog_tot_prior$Snipers, lwd = 1, lty = 6, col = 'gray60', add = TRUE)
dens(post_cog_tot_prior$`Military Controls`, lwd = 1, lty = 6, col = 'gray40', add = TRUE)
legend(-1.5, 3, box.lty = 0, cex = 0.9,
       legend = c("Breachers", "Snipers", "Military Controls"),
       fill = c("darkorange", "blue3", "forestgreen"))
text(-1.2, 0.37, "Breachers Prior", cex = 0.7, col = 'gray80')
text(-1.2, 0.5, "Snipers Prior", cex = 0.7, col = 'gray60')
text(-1.2, 0.17, "Military Controls Prior", cex = 0.7, col = 'gray40')
text(-1.5, 3, "A", cex = 0.9, font = 2)

## B. Unadjusted Contrasts (all on same plot)
plot(density(as.numeric(post_cog_tot_df$contrast_bc), adjust = 0.5),
     xlim = c(-1.5, 1.5), ylim = c(0, 4),
     lwd = 3, col = "firebrick", 
     main = "",
     font.lab = 2, cex.lab = 0.9, 
     xlab = "Difference in Latent Neurocognitive Scores (Unadjusted, SD Units)", 
     ylab = "", frame = FALSE)
lines(density(as.numeric(post_cog_tot_df$contrast_sc), adjust = 0.5),
      lwd = 3, col = "steelblue4")
lines(density(as.numeric(post_cog_tot_df$contrast_bs), adjust = 0.5),
      lwd = 3, col = "goldenrod3")
abline(v = 0, lty = 3, col = "gray40")
legend("topright", box.lty = 0, cex = 0.9,
       legend = c("Breachers - Military Controls", "Snipers - Military Controls", "Breachers - Snipers"),
       col = c("firebrick", "steelblue4", "goldenrod3"), lwd = 3)
text(-1.5, 4, "B", cex = 0.9, font = 2)

## C. Adjusted Group Estimates ptsd
dens(post_cog_dir_ptsd_df$Breachers, xlim = c(-1.5, 1.5),
     lwd = 3, col = "darkorange", 
     xlab = "Latent Neurocognitive Scores (PTSD Adjusted, SD Units)", 
     ylim = c(0, 3), font.lab = 2, 
     cex.lab = 0.9, frame = FALSE)
dens(post_cog_dir_ptsd_df$Snipers, lwd = 3, col = "blue3", add = TRUE)
dens(post_cog_dir_ptsd_df$`Military Controls`, lwd = 3, col = "forestgreen", add = TRUE)
dens(post_cog_dir_prior$Breachers, lwd = 1, lty = 6, col = 'gray80', add = TRUE)
dens(post_cog_dir_prior$Snipers, lwd = 1, lty = 6, col = 'gray60', add = TRUE)
dens(post_cog_dir_prior$`Military Controls`, lwd = 1, lty = 6, col = 'gray40', add = TRUE)
text(-1.5, 3, "C", cex = 0.9, font = 2)

## D. Adjusted Contrasts ptsd
plot(density(as.numeric(post_cog_dir_ptsd_df$contrast_bc), adjust = 0.5),
     xlim = c(-1.5, 1.5), ylim = c(0, 3), 
     lwd = 3, col = "firebrick",  
     font.lab = 2, cex.lab = 0.9,
     main = "",
     xlab = "Difference in Latent Neurocognitive Scores (PTSD Adjusted, SD Units)",
     ylab = "", frame = FALSE)
lines(density(as.numeric(post_cog_dir_ptsd_df$contrast_sc), adjust = 0.5),
      lwd = 3, col = "steelblue4")
lines(density(as.numeric(post_cog_dir_ptsd_df$contrast_bs), adjust = 0.5),
      lwd = 3, col = "goldenrod3")
abline(v = 0, lty = 3, col = "gray40")
text(-1.5, 3, "D", cex = 0.9, font = 2)

## E. Adjusted Group Estimates bsis
dens(post_cog_dir_bsis_df$Breachers, xlim = c(-1.5, 1.5),
     lwd = 3, col = "darkorange", 
     xlab = "Latent Neurocognitive Scores (BSI-S Adjusted, SD Units)", 
     ylim = c(0, 3), font.lab = 2, 
     cex.lab = 0.9, frame = FALSE)
dens(post_cog_dir_bsis_df$Snipers, lwd = 3, col = "blue3", add = TRUE)
dens(post_cog_dir_bsis_df$`Military Controls`, lwd = 3, col = "forestgreen", add = TRUE)
dens(post_cog_dir_prior$Breachers, lwd = 1, lty = 6, col = 'gray80', add = TRUE)
dens(post_cog_dir_prior$Snipers, lwd = 1, lty = 6, col = 'gray60', add = TRUE)
dens(post_cog_dir_prior$`Military Controls`, lwd = 1, lty = 6, col = 'gray40', add = TRUE)
text(-1.5, 3, "E", cex = 0.9, font = 2)

## F. Adjusted Contrasts bsis
plot(density(as.numeric(post_cog_dir_bsis_df$contrast_bc), adjust = 0.5),
     xlim = c(-1.5, 1.5), ylim = c(0, 3), 
     lwd = 3, col = "firebrick",  
     font.lab = 2, cex.lab = 0.9,
     main = "",
     xlab = "Difference in Latent Neurocognitive Scores (BSI-S Adjusted, SD Units)",
     ylab = "", frame = FALSE)
lines(density(as.numeric(post_cog_dir_bsis_df$contrast_sc), adjust = 0.5),
      lwd = 3, col = "steelblue4")
lines(density(as.numeric(post_cog_dir_bsis_df$contrast_bs), adjust = 0.5),
      lwd = 3, col = "goldenrod3")
abline(v = 0, lty = 3, col = "gray40")
text(-1.5, 3, "F", cex = 0.9, font = 2)

## G. Adjusted Group Estimates bsia
dens(post_cog_dir_bsia_df$Breachers, xlim = c(-1.5, 1.5),
     lwd = 3, col = "darkorange", 
     xlab = "Latent Neurocognitive Scores (BSI-A Adjusted, SD Units)", 
     ylim = c(0, 3), font.lab = 2, 
     cex.lab = 0.9, frame = FALSE)
dens(post_cog_dir_bsia_df$Snipers, lwd = 3, col = "blue3", add = TRUE)
dens(post_cog_dir_bsia_df$`Military Controls`, lwd = 3, col = "forestgreen", add = TRUE)
dens(post_cog_dir_prior$Breachers, lwd = 1, lty = 6, col = 'gray80', add = TRUE)
dens(post_cog_dir_prior$Snipers, lwd = 1, lty = 6, col = 'gray60', add = TRUE)
dens(post_cog_dir_prior$`Military Controls`, lwd = 1, lty = 6, col = 'gray40', add = TRUE)
text(-1.5, 3, "G", cex = 0.9, font = 2)

## H. Adjusted Contrasts bsia
plot(density(as.numeric(post_cog_dir_bsia_df$contrast_bc), adjust = 0.5),
     xlim = c(-1.5, 1.5), ylim = c(0, 3), 
     lwd = 3, col = "firebrick",  
     font.lab = 2, cex.lab = 0.9,
     main = "",
     xlab = "Difference in Latent Neurocognitive Scores (BSI-A Adjusted, SD Units)",
     ylab = "", frame = FALSE)
lines(density(as.numeric(post_cog_dir_bsia_df$contrast_sc), adjust = 0.5),
      lwd = 3, col = "steelblue4")
lines(density(as.numeric(post_cog_dir_bsia_df$contrast_bs), adjust = 0.5),
      lwd = 3, col = "goldenrod3")
abline(v = 0, lty = 3, col = "gray40")
text(-1.5, 3, "H", cex = 0.9, font = 2)

## I. Adjusted Group Estimates bsid
dens(post_cog_dir_bsid_df$Breachers, xlim = c(-1.5, 1.5),
     lwd = 3, col = "darkorange", 
     xlab = "Latent Neurocognitive Scores (BSI-D Adjusted, SD Units)", 
     ylim = c(0, 3), font.lab = 2, 
     cex.lab = 0.9, frame = FALSE)
dens(post_cog_dir_bsid_df$Snipers, lwd = 3, col = "blue3", add = TRUE)
dens(post_cog_dir_bsid_df$`Military Controls`, lwd = 3, col = "forestgreen", add = TRUE)
dens(post_cog_dir_prior$Breachers, lwd = 1, lty = 6, col = 'gray80', add = TRUE)
dens(post_cog_dir_prior$Snipers, lwd = 1, lty = 6, col = 'gray60', add = TRUE)
dens(post_cog_dir_prior$`Military Controls`, lwd = 1, lty = 6, col = 'gray40', add = TRUE)
text(-1.5, 3, "I", cex = 0.9, font = 2)

## J. Adjusted Contrasts bsid
plot(density(as.numeric(post_cog_dir_bsid_df$contrast_bc), adjust = 0.5),
     xlim = c(-1.5, 1.5), ylim = c(0, 3), 
     lwd = 3, col = "firebrick",  
     font.lab = 2, cex.lab = 0.9,
     main = "",
     xlab = "Difference in Latent Neurocognitive Scores (BSI-D Adjusted, SD Units)",
     ylab = "", frame = FALSE)
lines(density(as.numeric(post_cog_dir_bsid_df$contrast_sc), adjust = 0.5),
      lwd = 3, col = "steelblue4")
lines(density(as.numeric(post_cog_dir_bsid_df$contrast_bs), adjust = 0.5),
      lwd = 3, col = "goldenrod3")
abline(v = 0, lty = 3, col = "gray40")
text(-1.5, 3, "J", cex = 0.9, font = 2)

sup_fig1 <- recordPlot()
dev.off()

## PPC ##
group_vec <- as.factor(df$group)
group_levels <- levels(group_vec)
z_scores_raw <- t(apply(post_z_cog, 1, scale)) 

# Means by group across all posterior draws
empirical_means <- matrix(NA, nrow = nrow(z_scores_raw), ncol = length(group_levels))
colnames(empirical_means) <- group_levels

for (g in 1:length(group_levels)) {
  idx <- which(group_vec == group_levels[g])
  empirical_means[, g] <- rowMeans(z_scores_raw[, idx, drop = FALSE])
}

# Now calculate mean and 90% interval
empirical_summary <- data.frame(
  group = group_levels,
  mean = apply(empirical_means, 2, mean),
  lower = apply(empirical_means, 2, function(x) quantile(x, 0.05)),
  upper = apply(empirical_means, 2, function(x) quantile(x, 0.95))
)

# PPC Plot 
pdf("ppc_cog.pdf", width = 8, height = 6)
par(mfrow = c(3, 2), mar = c(4, 4, 3, 2))

## Panel A: Unadjusted
dens(post_cog_tot_df$Breachers, xlim = c(-1.5, 1.5), ylim = c(0, 3),
     col = "pink", lwd = 3, main = "",
     xlab = "Latent Neurocognitive Scores (Unadjusted, SD Units)", font.lab = 2, frame = FALSE)
dens(post_cog_tot_df$Snipers, add = TRUE, col = "purple", lwd = 3)
dens(post_cog_tot_df$`Military Controls`, add = TRUE, col = "lightgreen", lwd = 3)
abline(v = empirical_summary$mean[empirical_summary$group == "Breacher"],
       col = "pink", lty = 2, lwd = 2)
abline(v = empirical_summary$mean[empirical_summary$group == "Sniper"],
       col = "purple", lty = 2, lwd = 2)
abline(v = empirical_summary$mean[empirical_summary$group == "Military Control"],
       col = "lightgreen", lty = 2, lwd = 2)
legend("topright", bty = "n", cex = 0.9,
       legend = c("Posterior", "Empirical Mean"),
       lty = c(1, 2), lwd = 2, col = "black")
text(-1.5, 3, "A", cex = 0.9, font = 2)

## Panel B: Adjusted PTSD
dens(post_cog_dir_ptsd_df$Breachers, xlim = c(-1.5, 1.5), ylim = c(0, 4),
     col = "pink", lwd = 3, main = "",
     xlab = "Latent Neurocognitive Scores (PTSD Adjusted, SD Units)", font.lab = 2, frame = FALSE)
dens(post_cog_dir_ptsd_df$Snipers, add = TRUE, col = "purple", lwd = 3)
dens(post_cog_dir_ptsd_df$`Military Controls`, add = TRUE, col = "lightgreen", lwd = 3)
abline(v = empirical_summary$mean[empirical_summary$group == "Breacher"],
       col = "pink", lty = 2, lwd = 2)
abline(v = empirical_summary$mean[empirical_summary$group == "Sniper"],
       col = "purple", lty = 2, lwd = 2)
abline(v = empirical_summary$mean[empirical_summary$group == "Military Control"],
       col = "lightgreen", lty = 2, lwd = 2)
legend(0.6, 4.3, box.lty = 0, cex = 0.9,
       legend = c("Breachers", "Snipers", "Military Controls"),
       fill = c("pink", "purple", "lightgreen"))
text(-1.5, 4, "B", cex = 0.9, font = 2)

## Panel C: Adjusted bsis
dens(post_cog_dir_bsis_df$Breachers, xlim = c(-1.5, 1.5), ylim = c(0, 3),
     col = "pink", lwd = 3, main = "",
     xlab = "Latent Neurocognitive Scores (BSI-S Adjusted, SD Units)", font.lab = 2, frame = FALSE)
dens(post_cog_dir_bsis_df$Snipers, add = TRUE, col = "purple", lwd = 3)
dens(post_cog_dir_bsis_df$`Military Controls`, add = TRUE, col = "lightgreen", lwd = 3)
abline(v = empirical_summary$mean[empirical_summary$group == "Breacher"],
       col = "pink", lty = 2, lwd = 2)
abline(v = empirical_summary$mean[empirical_summary$group == "Sniper"],
       col = "purple", lty = 2, lwd = 2)
abline(v = empirical_summary$mean[empirical_summary$group == "Military Control"],
       col = "lightgreen", lty = 2, lwd = 2)
text(-1.5, 3, "C", cex = 0.9, font = 2)

## Panel D: Adjusted bsia
dens(post_cog_dir_bsia_df$Breachers, xlim = c(-1.5, 1.5), ylim = c(0, 3),
     col = "pink", lwd = 3, main = "",
     xlab = "Latent Neurocognitive Scores (BSI-A Adjusted, SD Units)", font.lab = 2, frame = FALSE)
dens(post_cog_dir_bsia_df$Snipers, add = TRUE, col = "purple", lwd = 3)
dens(post_cog_dir_bsia_df$`Military Controls`, add = TRUE, col = "lightgreen", lwd = 3)
abline(v = empirical_summary$mean[empirical_summary$group == "Breacher"],
       col = "pink", lty = 2, lwd = 2)
abline(v = empirical_summary$mean[empirical_summary$group == "Sniper"],
       col = "purple", lty = 2, lwd = 2)
abline(v = empirical_summary$mean[empirical_summary$group == "Military Control"],
       col = "lightgreen", lty = 2, lwd = 2)
text(-1.5, 3, "D", cex = 0.9, font = 2)

## Panel E: Adjusted bsid
dens(post_cog_dir_bsid_df$Breachers, xlim = c(-1.5, 1.5), ylim = c(0, 3),
     col = "pink", lwd = 3, main = "",
     xlab = "Latent Neurocognitive Scores (BSI-D Adjusted, SD Units)", font.lab = 2, frame = FALSE)
dens(post_cog_dir_bsid_df$Snipers, add = TRUE, col = "purple", lwd = 3)
dens(post_cog_dir_bsid_df$`Military Controls`, add = TRUE, col = "lightgreen", lwd = 3)
abline(v = empirical_summary$mean[empirical_summary$group == "Breacher"],
       col = "pink", lty = 2, lwd = 2)
abline(v = empirical_summary$mean[empirical_summary$group == "Sniper"],
       col = "purple", lty = 2, lwd = 2)
abline(v = empirical_summary$mean[empirical_summary$group == "Military Control"],
       col = "lightgreen", lty = 2, lwd = 2)
text(-1.5, 3, "E", cex = 0.9, font = 2)

ppc_cog_df <- recordPlot()
dev.off()

#for results reporting
precis(post_cog_tot_df, 2, probs = 0.9)
mean(post_cog_tot_df$contrast_bc >0)
mean(post_cog_tot_df$contrast_sc >0)
mean(post_cog_tot_df$contrast_bs >0)
precis(post_cog_dir_ptsd_df, 2, probs = 0.9)
precis(post_cog_dir_bsis_df, 2, probs = 0.9)
precis(post_cog_dir_bsia_df, 2, probs = 0.9)
precis(post_cog_dir_bsid_df, 2, probs = 0.9)

# Raw data sanity checks
aggregate(colMeans(post_z_cog), by = list(group = df$group), FUN = mean)
# Transpose so rows = subjects, cols = posterior draws
post_z_subject <- t(post_z_cog)

# Create vector of group labels
group_vec <- df$group

# Compute mean posterior latent score for each subject
subject_means <- rowMeans(post_z_subject)

# Now group-level means:
tapply(subject_means, group_vec, mean)


