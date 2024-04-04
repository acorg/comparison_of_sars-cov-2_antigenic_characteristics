
rm(list = ls())

library(Racmacs)
library(cmdstanr)
library(titertools)
library(bayesplot)
library(tidyverse)

source('code/data_generation/load_maps_for_comparison.R')
source('code/functions/gmts.R')

# Load the maps
maps <- load_maps_for_comparison('', sr_group_standardisation = 'data/metadata/sr_group_standardisation_WT_conv_as_D614G_as614D.csv')

# Re-arrange the data
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
  sr_individual_effects = rep(0, numSera(data)),
  map_effects = rep(0, numLayers(data)),
  logtiter_error_sigma = rep(0.8, numLayers(data))
)

# Fetch the model
mod <- cmdstan_model('code/stan_models/titer_magnitude.stan')

# Run the model
mod_optimized <- mod$sample(
  data = data_list,
  init = list(data_init, data_init, data_init, data_init),
  chains = 4, 
  refresh = 500,
  iter_sampling  = 10000,
  iter_warmup = 1000
)

# Draw from the posterior
draws <- mod_optimized$draws(format = 'df', variables = c('dataset_effects'))

# Re-name columns
colnames(draws) <- c('duke', 'maryland', 'galveston', 'emory',
                     'madison_pooled', 'madison_unpooled', 'st_louis',
                     'oxford', 'mt_sinai_human', 'emc_prnt',
                     'innsbruck', 'charite', 'madison_frnt', 'fda', 'geneva',
                     'amc', 'emc_calu', 'emc_vero', '.chain', '.iteration', '.draw')


# Save the draws
saveRDS(draws, 'data/titer_analyses/titer_magnitude/dataset_magnitude_effect/dataset_magnitude_effect_posterior_samples.rds')


# Check convergence
mod_optimized$summary('dataset_effects')
# # A tibble: 18 × 10
# variable               mean  median    sd   mad    q5     q95  rhat ess_bulk ess_tail
# <chr>                 <dbl>   <dbl> <dbl> <dbl> <dbl>   <dbl> <dbl>    <dbl>    <dbl>
#   1 dataset_effects[1]  -3.71   -3.71    1.05  1.04 -5.41 -2.00    1.00    1291.    2629.
# 2 dataset_effects[2]   2.56    2.56    2.57  2.59 -1.68  6.80    1.00    8940.   18252.
# 3 dataset_effects[3]  -0.287  -0.268   1.73  1.72 -3.14  2.52    1.00    3198.    6228.
# 4 dataset_effects[4]  -3.93   -3.93    1.22  1.21 -5.96 -1.90    1.00    1668.    2974.
# 5 dataset_effects[5]  -2.25   -2.26    2.15  2.16 -5.80  1.31    1.00    6079.   13294.
# 6 dataset_effects[6]  -2.42   -2.42    1.53  1.53 -4.95  0.0624  1.00    2617.    6844.
# 7 dataset_effects[7]   0.858   0.868   1.37  1.38 -1.42  3.11    1.00    1869.    3862.
# 8 dataset_effects[8]  -3.41   -3.40    1.07  1.07 -5.19 -1.66    1.00    1184.    2101.
# 9 dataset_effects[9]  -4.03   -4.04    1.46  1.45 -6.44 -1.63    1.00    1935.    4334.
# 10 dataset_effects[10] -0.224  -0.234   1.51  1.51 -2.68  2.26    1.00    2903.    6494.
# 11 dataset_effects[11] -3.65   -3.63    1.20  1.22 -5.66 -1.73    1.00    1635.    3024.
# 12 dataset_effects[12] -0.589  -0.588   1.61  1.61 -3.27  2.08    1.00    3131.    6648.
# 13 dataset_effects[13]  0.0337  0.0428  1.51  1.51 -2.44  2.52    1.00    2824.    5543.
# 14 dataset_effects[14] -2.38   -2.37    1.15  1.16 -4.28 -0.498   1.00    1530.    3559.
# 15 dataset_effects[15] -6.69   -6.69    1.11  1.11 -8.51 -4.89    1.00    1560.    3694.
# 16 dataset_effects[16] -4.57   -4.58    1.31  1.31 -6.74 -2.45    1.00    1785.    3668.
# 17 dataset_effects[17] -0.595  -0.591   1.60  1.58 -3.23  2.01    1.00    2901.    6203.
# 18 dataset_effects[18] -0.934  -0.939   1.59  1.58 -3.54  1.68    1.00    2724.    6217.

# Look at traces
mcmc_trace(mod_optimized$draws(variables = c('dataset_effects'))) -> gp
ggsave('data/titer_analyses/titer_magnitude/dataset_magnitude_effect/dataset_magnitude_effect_posterior_samples_traces.png', width = 15, height = 12)
