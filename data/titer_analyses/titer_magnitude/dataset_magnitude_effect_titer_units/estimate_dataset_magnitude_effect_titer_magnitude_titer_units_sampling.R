
rm(list = ls())

library(Racmacs)
library(cmdstanr)
library(bayesplot)
library(titertools)
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
  N_titer_units = 4,
  # NT50: 1, NT90: 2, NT99: 3, IU/ml: 4
  titer_units = c(1, 3, 1, 1, 3, 3, 1, 1, 1, 1, 1, 2, 1, 1, 2, 4, 1, 1),
  logtiters = simplify2array(logtiter_layers),
  titertypes = simplify2array(Racmacs:::titertypesTableLayers(data))
)

# Decide on starting conditions for the model
sr_group_gmts <- srGroupGMTs(data, ci_method = 'quap')
sr_group_gmts[is.na(sr_group_gmts)] <- 7

data_init <- list(
  sr_group_gmts = sr_group_gmts,
  sr_effects = rep(0, numSera(data)),
  animal_effects = rep(0, numLayers(data)),
  assay_effects = rep(0, numLayers(data)),
  logtiter_error_sigma = rep(0.8, numLayers(data))
)

# Fetch the model
mod <- cmdstan_model('code/stan_models/dataset_magnitude_titer_magnitude_titer_units.stan')

# Run the model
mod_optimized <- mod$sample(
  data = data_list,
  init = list(data_init, data_init, data_init, data_init),
  chains = 4, 
  refresh = 500,
  iter_sampling  = 10000,
  iter_warmup = 1000
)

# Warning: 39939 of 40000 (100.0%) transitions hit the maximum treedepth limit of 10.

# Draw from the posterior
draws_titer_units <- mod_optimized$draws(format = 'df', variables = c('titer_units_effects'))
draws_dataset_effects <- mod_optimized$draws(format = 'df', variables = c('dataset_effects'))

# Re-name columns
colnames(draws_titer_units) <- c('NT50', 'NT90', 'NT99', 'IU/ml', '.chain', '.iteration', '.draw')
colnames(draws_dataset_effects) <- c(
  'duke', 'maryland', 'galveston', 'emory', 'madison_pooled', 'madison_unpooled',
  'st_louis', 'oxford', 'mt_sinai_human', 'emc_prnt', 'innsbruck', 'charite',
  'madison_frnt', 'fda', 'geneva', 'amc', 'emc_calu', 'emc_vero',
  '.chain', '.iteration', '.draw')


# Save the draws
saveRDS(draws_titer_units, 'data/titer_analyses/titer_magnitude/dataset_magnitude_effect_titer_units/dataset_magnitude_effect_titer_magnitude_titer_units_titer_units_draws.rds')
saveRDS(draws_dataset_effects, 'data/titer_analyses/titer_magnitude/dataset_magnitude_effect_titer_units/dataset_magnitude_effect_titer_magnitude_titer_units_dataset_magnitude_draws.rds')


# Check convergence
mod_optimized$summary('titer_units_effects')
# # A tibble: 4 × 10
# variable                mean median    sd   mad    q5   q95  rhat ess_bulk ess_tail
# <chr>                  <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>    <dbl>    <dbl>
#   1 titer_units_effects[1] -2.69  -2.70  1.96  1.94 -5.89 0.539  1.00    2130.    6023.
# 2 titer_units_effects[2] -3.11  -3.11  3.59  3.58 -8.99 2.83   1.00   13041.   22752.
# 3 titer_units_effects[3] -1.26  -1.24  3.19  3.17 -6.56 3.99   1.00    9027.   23384.
# 4 titer_units_effects[4] -2.82  -2.83  4.30  4.30 -9.88 4.26   1.00   29021.   26758.

mod_optimized$summary('dataset_effects')
# # A tibble: 18 × 10
# variable              mean median    sd   mad      q5   q95  rhat ess_bulk ess_tail
# <chr>                <dbl>  <dbl> <dbl> <dbl>   <dbl> <dbl> <dbl>    <dbl>    <dbl>
# 1 dataset_effects[1]  -1.99  -2.00   1.77  1.77  -4.86  0.944  1.00    7410.   14571.
# 2 dataset_effects[2]   2.79   2.81   3.57  3.56  -3.06  8.66   1.00   31120.   30254.
# 3 dataset_effects[3]   1.37   1.35   2.21  2.20  -2.27  5.01   1.00    9624.   16438.
# 4 dataset_effects[4]  -2.28  -2.30   1.87  1.85  -5.32  0.796  1.00    7199.   14573.
# 5 dataset_effects[5]  -1.95  -1.94   3.42  3.42  -7.55  3.69   1.00   33895.   30828.
# 6 dataset_effects[6]  -2.10  -2.11   3.22  3.22  -7.39  3.25   1.00   22960.   27914.
# 7 dataset_effects[7]   2.48   2.47   1.95  1.94  -0.725 5.70   1.00    6923.   12909.
# 8 dataset_effects[8]  -1.70  -1.70   1.79  1.78  -4.64  1.26   1.00    6847.   14524.
# 9 dataset_effects[9]  -2.38  -2.39   1.98  1.99  -5.62  0.898  1.00    5982.   12528.
# 10 dataset_effects[10]  1.37   1.35   2.05  2.06  -2.01  4.74   1.00    8603.   15771.
# 11 dataset_effects[11] -2.00  -2.02   1.85  1.85  -5.00  1.05   1.00    7384.   14936.
# 12 dataset_effects[12]  1.44   1.44   3.61  3.60  -4.51  7.44   1.00   27850.   28001.
# 13 dataset_effects[13]  1.64   1.65   2.05  2.03  -1.75  5.02   1.00    8274.   14812.
# 14 dataset_effects[14] -0.689 -0.694  1.84  1.83  -3.70  2.34   1.00    7482.   14987.
# 15 dataset_effects[15] -4.59  -4.58   3.56  3.56 -10.4   1.31   1.00   25282.   27353.
# 16 dataset_effects[16] -2.80  -2.79   4.29  4.28  -9.84  4.24   1.00   30577.   27449.
# 17 dataset_effects[17]  1.01   1.01   2.07  2.08  -2.40  4.38   1.00    8975.   15103.
# 18 dataset_effects[18]  0.658  0.662  2.07  2.09  -2.77  4.04   1.00    8510.   15943.

# Look at traces
mcmc_trace(mod_optimized$draws(variables = c('titer_units_effects'))) -> gp
ggsave('data/titer_analyses/titer_magnitude/dataset_magnitude_effect_titer_units/traces-dataset-effects-titer-units.png', width = 15, height = 12)

mcmc_trace(mod_optimized$draws(variables = c('dataset_effects'))) -> gp
ggsave('data/titer_analyses/titer_magnitude/dataset_magnitude_effect_titer_units/traces-dataset-effects.png', width = 15, height = 12)
