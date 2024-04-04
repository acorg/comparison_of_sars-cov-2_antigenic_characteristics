
rm(list = ls())

library(Racmacs)
library(cmdstanr)
library(bayesplot)
library(titertools)

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
  N_animals = 3,
  N_assays = 6,
  # human: 1, hamster: 2, mouse: 3
  animals = c(1, 3, 2, 1, 2, 2, 3, 1, 1, 2, 1, 2, 2, 1, 1, 1, 2, 2),
  # FRNT: 1, PV-neut: 2, VSV-PV-neut: 3, PRNT: 4, Microneut: 5, CPE: 6
  assays = c(2, 6, 1, 1, 6, 6, 1, 1, 5, 4, 1, 4, 1, 2, 4, 2, 3, 3),
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
mod <- cmdstan_model('code/stan_models/titer_magnitude_animal_model_assay_differences.stan')

# Run the model
mod_optimized <- mod$sample(
  data = data_list,
  init = list(data_init, data_init, data_init, data_init),
  chains = 4, 
  refresh = 500,
  iter_sampling  = 10000,
  iter_warmup = 1000
)

# Warning: 8286 of 40000 (21.0%) transitions hit the maximum treedepth limit of 10.

# Draw from the posterior
draws_animal <- mod_optimized$draws(format = 'df', variables = c('animal_effects'))
draws_assay <- mod_optimized$draws(format = 'df', variables = c('assay_effects'))

# Re-name columns
colnames(draws_animal) <- c('human', 'hamster', 'mouse', '.chain', '.iteration', '.draw')
colnames(draws_assay) <- c('FRNT', 'PV-neut', 'VSV-PV-neut', 'PRNT', 'Microneut', 'CPE', '.chain', '.iteration', '.draw')

# Save the draws
saveRDS(draws_animal, 'data/titer_analyses/titer_magnitude/dataset_magnitude_effect_animal_model_assay_differences/dataset_magnitude_effect_animal_model_differences_draws.rds')
saveRDS(draws_assay, 'data/titer_analyses/titer_magnitude/dataset_magnitude_effect_animal_model_assay_differences/dataset_magnitude_effect_assay_differences_draws.rds')

# Check convergence
mod_optimized$summary('animal_effects')
# # A tibble: 3 × 10
# variable            mean median    sd   mad    q5   q95  rhat ess_bulk ess_tail
# <chr>              <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>    <dbl>    <dbl>
# 1 animal_effects[1] -5.16  -5.16   2.14  2.16 -8.69 -1.61  1.00    1067.    4860.
# 2 animal_effects[2] -0.754 -0.763  2.16  2.18 -4.31  2.83  1.00    1437.    4563.
# 3 animal_effects[3]  0.177  0.169  2.26  2.30 -3.51  3.91  1.00    1063.    4437.

mod_optimized$summary('assay_effects')
# # A tibble: 6 × 10
# variable             mean  median    sd   mad    q5   q95  rhat ess_bulk ess_tail
# <chr>               <dbl>   <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>    <dbl>    <dbl>
# 1 assay_effects[1]  0.0395  0.0398  2.06  2.06 -3.38  3.40  1.00    2200.    6394.
# 2 assay_effects[2]  0.388   0.373   2.09  2.09 -3.02  3.82  1.00    1919.    5479.
# 3 assay_effects[3] -1.28   -1.28    2.23  2.22 -4.95  2.37  1.00    1661.    4413.
# 4 assay_effects[4] -2.19   -2.19    2.09  2.09 -5.66  1.26  1.00    2283.    5709.
# 5 assay_effects[5] -0.130  -0.132   2.29  2.29 -3.92  3.63  1.00    1349.    4169.
# 6 assay_effects[6] -2.21   -2.21    2.24  2.23 -5.91  1.48  1.00    2270.    5182.

# Look at traces
mcmc_trace(mod_optimized$draws(variables = c('animal_effects', 'assay_effects'))) -> gp
ggsave('data/titer_analyses/titer_magnitude/dataset_magnitude_effect_animal_model_assay_differences/traces.png', width = 15, height = 12)
