
rm(list = ls())

library(tidyverse)
library(Racmacs)
library(cmdstanr)
library(titertools)
library(bayesplot)

source('code/data_generation/load_maps_for_comparison.R')
source('code/functions/gmts.R')

# Read in data
maps <- load_maps_for_comparison('', sr_group_standardisation = 'data/metadata/sr_group_standardisation_WT_conv_as_D614G_as614D.csv')

# Arrange the data
data <- arrange_data(maps)

# Run the model
options(mc.cores = parallel::detectCores())

# Replace logtiter table layer NAs with 0s, they will be dealt with by titer types
logtiter_layers <- Racmacs:::logtiterTableLayers(data)
logtiter_layers <- lapply(logtiter_layers, \(logtiters) {
  logtiters[is.na(logtiters)] <- 0
  logtiters
})

# Assemble the data list for the model
data_list <- list(
  only_prior = 0,
  N_datasets = numLayers(data),
  N_ags = numAntigens(data),
  N_srs = numSera(data),
  N_sr_groups = numSeraGroups(data),
  sr_groups = as.numeric(srGroups(data)),
  logtiters = simplify2array(logtiter_layers),
  titertypes = simplify2array(Racmacs:::titertypesTableLayers(data))
)

# Decide on starting conditions for the model
sr_group_gmts <- srGroupGMTs(data, ci_method = 'quap')
sr_group_gmts[is.na(sr_group_gmts)] <- 7

data_init <- list(
  sr_group_gmts = sr_group_gmts,
  sr_effects = rep(0, numSera(data)),
  map_effects = rep(0, numLayers(data)),
  logtiter_error_sigma = rep(0.8, numLayers(data))
)

# Fetch the model
mod <- cmdstan_model(stan_file = 'code/stan_models/titer_variability.stan')

mod_optimized <- mod$sample(
  data = data_list,
  init = list(data_init, data_init, data_init, data_init),
  chains = 4, 
  refresh = 500,
  iter_sampling  = 10000,
  iter_warmup = 1000
)

# Warning: 59 of 40000 (0.0%) transitions ended with a divergence.
# See https://mc-stan.org/misc/warnings for details.
# 
# Warning: 4 of 4 chains had an E-BFMI less than 0.2.
# See https://mc-stan.org/misc/warnings for details.

# Draw from the posterior
mod_optimized_draws <- mod_optimized$draws(format = 'df', variables = c('sr_group_gmts_deviation_sd'))

# Re-name columns
colnames(mod_optimized_draws) <- c('duke', 'maryland', 'galveston', 'emory',
                                   'madison_pooled', 'madison_unpooled', 'st_louis',
                                   'oxford', 'mt_sinai_human', 'emc_prnt',
                                   'innsbruck', 'charite', 'madison_frnt', 'fda', 'geneva',
                                   'amc', 'emc_calu', 'emc_vero', '.chain', '.iteration', '.draw')

# Save the draws
saveRDS(mod_optimized_draws, 'data/titer_analyses/titer_variability/titer_variability_effect/sr_group_gmts_deviation_sd_draws.rds')


# Check convergence
mod_optimized$summary('sr_group_gmts_deviation_sd')
# # A tibble: 18 × 10
# variable                        mean median     sd    mad    q5   q95  rhat ess_bulk ess_tail
# <chr>                          <dbl>  <dbl>  <dbl>  <dbl> <dbl> <dbl> <dbl>    <dbl>    <dbl>
#   1 sr_group_gmts_deviation_sd[1]  0.999  0.995 0.0991 0.0985 0.844 1.17   1.00    3127.    6451.
# 2 sr_group_gmts_deviation_sd[2]  1.78   1.79  0.299  0.267  1.28  2.26   1.01     736.     534.
# 3 sr_group_gmts_deviation_sd[3]  0.597  0.579 0.144  0.138  0.396 0.867  1.01     471.    1050.
# 4 sr_group_gmts_deviation_sd[4]  0.589  0.583 0.0888 0.0877 0.454 0.744  1.00    1219.    2207.
# 5 sr_group_gmts_deviation_sd[5]  0.964  0.963 0.179  0.177  0.675 1.26   1.01     745.     863.
# 6 sr_group_gmts_deviation_sd[6]  0.873  0.862 0.141  0.139  0.660 1.12   1.00    1166.    2464.
# 7 sr_group_gmts_deviation_sd[7]  1.30   1.23  0.447  0.391  0.727 2.13   1.01     239.     493.
# 8 sr_group_gmts_deviation_sd[8]  0.885  0.875 0.132  0.131  0.685 1.12   1.00    1375.    2681.
# 9 sr_group_gmts_deviation_sd[9]  0.849  0.815 0.221  0.205  0.556 1.26   1.01     487.     900.
# 10 sr_group_gmts_deviation_sd[10] 1.06   1.05  0.143  0.138  0.839 1.30   1.00    1540.    3529.
# 11 sr_group_gmts_deviation_sd[11] 0.866  0.857 0.122  0.117  0.683 1.08   1.00    1300.    3066.
# 12 sr_group_gmts_deviation_sd[12] 0.687  0.677 0.138  0.134  0.481 0.930  1.00     597.    1168.
# 13 sr_group_gmts_deviation_sd[13] 0.910  0.887 0.229  0.215  0.589 1.31   1.02     472.     770.
# 14 sr_group_gmts_deviation_sd[14] 0.851  0.846 0.104  0.103  0.690 1.03   1.00    1995.    4200.
# 15 sr_group_gmts_deviation_sd[15] 0.990  0.978 0.148  0.143  0.771 1.26   1.01    1256.    2173.
# 16 sr_group_gmts_deviation_sd[16] 0.738  0.725 0.144  0.139  0.528 0.994  1.01     664.    1601.
# 17 sr_group_gmts_deviation_sd[17] 0.601  0.589 0.118  0.115  0.428 0.813  1.01     783.    1760.
# 18 sr_group_gmts_deviation_sd[18] 0.861  0.850 0.167  0.163  0.608 1.15   1.00     766.    1397.

# Look at traces
mcmc_trace(mod_optimized$draws(variables = c('sr_group_gmts_deviation_sd'))) -> gp
ggsave('data/titer_analyses/titer_variability/titer_variability_effect/sr_group_gmts_deviation_sd_draws_traces.png', width = 15, height = 12)

