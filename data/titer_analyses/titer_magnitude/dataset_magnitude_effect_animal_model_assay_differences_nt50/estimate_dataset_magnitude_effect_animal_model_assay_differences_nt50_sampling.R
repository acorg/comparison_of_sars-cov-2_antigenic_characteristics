
rm(list = ls())

library(Racmacs)
library(cmdstanr)
library(bayesplot)
library(titertools)

source('code/data_generation/load_maps_for_comparison.R')
source('code/functions/gmts.R')

# Load the maps
maps <- load_maps_for_comparison('', sr_group_standardisation = 'data/metadata/sr_group_standardisation_WT_conv_as_D614G_as614D.csv')

# Remove the non-NT50 maps
maps <- maps[-c(2, 5, 6, 12, 15, 16)]

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
  animals = c(1, 2, 1, 3, 1, 1, 2, 1, 2, 1, 2, 2),
  # FRNT: 1, PV-neut: 2, VSV-PV-neut: 3, PRNT: 4, Microneut: 5, CPE: 6
  assays = c(2, 1, 1, 1, 1, 5, 4, 1, 1, 2, 3, 3),
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


# Draw from the posterior
draws_animal <- mod_optimized$draws(format = 'df', variables = c('animal_effects'))
draws_assay <- mod_optimized$draws(format = 'df', variables = c('assay_effects'))

# Re-name columns
colnames(draws_animal) <- c('human', 'hamster', 'mouse', '.chain', '.iteration', '.draw')
colnames(draws_assay) <- c('FRNT', 'PV-neut', 'VSV-PV-neut', 'PRNT', 'Microneut', 'CPE', '.chain', '.iteration', '.draw')

# Save the draws
saveRDS(draws_animal, 'data/titer_analyses/titer_magnitude/dataset_magnitude_effect_animal_model_assay_differences_nt50/estimate_dataset_magnitude_effect_animal_model_assay_differences_nt5_draws_animal.rds')
saveRDS(draws_assay, 'data/titer_analyses/titer_magnitude/dataset_magnitude_effect_animal_model_assay_differences_nt50/estimate_dataset_magnitude_effect_animal_model_assay_differences_nt5_draws_assay.rds')

# Check convergence
mod_optimized$summary('animal_effects')
# # A tibble: 3 × 10
# variable            mean median    sd   mad    q5   q95  rhat ess_bulk ess_tail
# <chr>              <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>    <dbl>    <dbl>
# 1 animal_effects[1] -4.18  -4.18   2.30  2.31 -7.95 -0.400  1.00    1710.    6648.
# 2 animal_effects[2] -0.790 -0.799  2.35  2.35 -4.64  3.08   1.00    1702.    6550.
# 3 animal_effects[3]  0.438  0.425  2.47  2.48 -3.60  4.53   1.00    1462.    5069.

mod_optimized$summary('assay_effects')
# # A tibble: 6 × 10
# variable             mean  median    sd   mad    q5   q95  rhat ess_bulk ess_tail
# <chr>               <dbl>   <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>    <dbl>    <dbl>
# 1 assay_effects[1] -0.827  -0.841   2.23  2.24 -4.47  2.84  1.00    2962.    8829.
# 2 assay_effects[2] -0.452  -0.462   2.28  2.31 -4.17  3.30  1.00    2446.    8643.
# 3 assay_effects[3] -1.46   -1.48    2.43  2.41 -5.46  2.56  1.00    2412.    6895.
# 4 assay_effects[4] -0.901  -0.914   2.53  2.50 -5.05  3.28  1.00    3565.    7984.
# 5 assay_effects[5] -0.956  -0.949   2.45  2.44 -4.98  3.06  1.00    1995.    5842.
# 6 assay_effects[6] -0.0634 -0.0792  5.96  5.91 -9.82  9.82  1.00   35955.   27889.

# Look at traces
mcmc_trace(mod_optimized$draws(variables = c('animal_effects', 'assay_effects'))) -> gp
ggsave('data/titer_analyses/titer_magnitude/dataset_magnitude_effect_animal_model_assay_differences_nt50/estimate_dataset_magnitude_effect_animal_model_assay_differences_nt5_traces.png', width = 15, height = 12)
