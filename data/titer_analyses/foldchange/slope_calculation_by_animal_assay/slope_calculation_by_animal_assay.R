
rm(list = ls())

library(tidyverse)
library(Racmacs)
library(cmdstanr)
library(patchwork)
library(titertools)
library(bayesplot)

source('code/data_generation/load_maps_for_comparison.R')
source('code/functions/map_longinfo.R')
source('code/metadata/common.R')

# Homologous ag info
homologous_ags <- c(
  'mRNA-1273 duke' = 'D614G',
  'mRNA-1273 maryland' = '614D',
  'mRNA-1273 galveston' = '614D',
  'mRNA-1273 emory' = 'D614G',
  'mRNA-1273 madison_pooled' = 'D614G',
  'mRNA-1273 madison_unpooled' = 'D614G',
  'mRNA-1273 st_louis' = 'D614G',
  'mRNA-1273 oxford' = '614D',
  'mRNA-1273 mt_sinai_mouse' = '614D',
  'mRNA-1273 mt_sinai_human' = '614D',
  'mRNA-1273 emc_prnt' = 'D614G',
  'mRNA-1273 kcl' = '614D',
  'mRNA-1273 innsbruck' = 'D614G',
  'mRNA-1273 fda' = 'D614G',
  'mRNA-1273 geneva' = 'D614G',
  'mRNA-1273 amc' = 'D614G',
  'mRNA-1273 charite' = 'D614G',
  'mRNA-1273 madison_frnt' = '614D',
  'mRNA-1273 emc_calu' = 'D614G',
  'mRNA-1273 emc_vero' = 'D614G',
  'D614G convalescent duke' = 'D614G',
  'D614G convalescent maryland' = '614D',
  'D614G convalescent galveston' = '614D',
  'D614G convalescent emory' = 'D614G',
  'D614G convalescent madison_pooled' = 'D614G',
  'D614G convalescent madison_unpooled' = 'D614G',
  'D614G convalescent st_louis' = 'D614G',
  'D614G convalescent oxford' = '614D',
  'D614G convalescent mt_sinai_mouse' = '614D',
  'D614G convalescent mt_sinai_human' = '614D',
  'D614G convalescent emc_prnt' = 'D614G',
  'D614G convalescent kcl' = '614D',
  'D614G convalescent innsbruck' = 'D614G',
  'D614G convalescent fda' = 'D614G',
  'D614G convalescent geneva' = 'D614G',
  'D614G convalescent amc' = 'D614G',
  'D614G convalescent charite' = 'D614G',
  'D614G convalescent madison_frnt' = '614D',
  'D614G convalescent emc_calu' = 'D614G',
  'D614G convalescent emc_vero' = 'D614G',
  'B.1.351 convalescent duke' = 'B.1.351',
  'B.1.351 convalescent maryland' = 'B.1.351',
  'B.1.351 convalescent galveston' = 'B.1.351',
  'B.1.351 convalescent emory' = 'B.1.351',
  'B.1.351 convalescent madison_pooled' = 'B.1.351',
  'B.1.351 convalescent madison_unpooled' = 'B.1.351',
  'B.1.351 convalescent st_louis' = 'B.1.351',
  'B.1.351 convalescent oxford' = 'B.1.351',
  'B.1.351 convalescent mt_sinai_mouse' = 'B.1.351',
  'B.1.351 convalescent mt_sinai_human' = 'B.1.351',
  'B.1.351 convalescent emc_prnt' = 'B.1.351',
  'B.1.351 convalescent kcl' = 'B.1.351',
  'B.1.351 convalescent innsbruck' = 'B.1.351',
  'B.1.351 convalescent fda' = 'B.1.351',
  'B.1.351 convalescent geneva' = 'B.1.351',
  'B.1.351 convalescent amc' = 'B.1.351',
  'B.1.351 convalescent charite' = 'B.1.351',
  'B.1.351 convalescent madison_frnt' = 'B.1.351',
  'B.1.351 convalescent emc_calu' = 'B.1.351',
  'B.1.351 convalescent emc_vero' = 'B.1.351',
  'B.1.617.2 convalescent duke' = 'B.1.617.2',
  'B.1.617.2 convalescent maryland' = 'B.1.617.2',
  'B.1.617.2 convalescent galveston' = 'B.1.617.2',
  'B.1.617.2 convalescent emory' = 'B.1.617.2',
  'B.1.617.2 convalescent madison_pooled' = 'B.1.617.2',
  'B.1.617.2 convalescent madison_unpooled' = 'B.1.617.2',
  'B.1.617.2 convalescent st_louis' = 'B.1.617.2',
  'B.1.617.2 convalescent oxford' = 'B.1.617.2',
  'B.1.617.2 convalescent mt_sinai_mouse' = 'B.1.617.2',
  'B.1.617.2 convalescent mt_sinai_human' = 'B.1.617.2',
  'B.1.617.2 convalescent emc_prnt' = 'B.1.617.2',
  'B.1.617.2 convalescent kcl' = 'B.1.617.2',
  'B.1.617.2 convalescent innsbruck' = 'B.1.617.2',
  'B.1.617.2 convalescent fda' = 'B.1.617.2',
  'B.1.617.2 convalescent geneva' = 'B.1.617.2',
  'B.1.617.2 convalescent amc' = 'B.1.617.2',
  'B.1.617.2 convalescent charite' = 'B.1.617.2',
  'B.1.617.2 convalescent madison_frnt' = 'B.1.617.2',
  'B.1.617.2 convalescent emc_calu' = 'B.1.617.2',
  'B.1.617.2 convalescent emc_vero' = 'B.1.617.2',
  'P.1 convalescent duke' = 'P.1',
  'P.1 convalescent maryland' = 'P.1',
  'P.1 convalescent galveston' = 'P.1',
  'P.1 convalescent emory' = 'P.1',
  'P.1 convalescent madison_pooled' = 'P.1',
  'P.1 convalescent madison_unpooled' = 'P.1',
  'P.1 convalescent st_louis' = 'P.1',
  'P.1 convalescent oxford' = 'P.1',
  'P.1 convalescent mt_sinai_mouse' = 'P.1',
  'P.1 convalescent mt_sinai_human' = 'P.1',
  'P.1 convalescent emc_prnt' = 'P.1',
  'P.1 convalescent kcl' = 'P.1',
  'P.1 convalescent innsbruck' = 'P.1',
  'P.1 convalescent fda' = 'P.1',
  'P.1 convalescent geneva' = 'P.1',
  'P.1 convalescent amc' = 'P.1',
  'P.1 convalescent charite' = 'P.1',
  'P.1 convalescent madison_frnt' = 'P.1',
  'P.1 convalescent emc_calu' = 'P.1',
  'P.1 convalescent emc_vero' = 'P.1',
  'B.1.1.7 convalescent duke' = 'B.1.1.7',
  'B.1.1.7 convalescent maryland' = 'B.1.1.7',
  'B.1.1.7 convalescent galveston' = 'B.1.1.7',
  'B.1.1.7 convalescent emory' = 'B.1.1.7',
  'B.1.1.7 convalescent madison_pooled' = 'B.1.1.7',
  'B.1.1.7 convalescent madison_unpooled' = 'B.1.1.7',
  'B.1.1.7 convalescent st_louis' = 'B.1.1.7',
  'B.1.1.7 convalescent oxford' = 'B.1.1.7',
  'B.1.1.7 convalescent mt_sinai_mouse' = 'B.1.1.7',
  'B.1.1.7 convalescent mt_sinai_human' = 'B.1.1.7',
  'B.1.1.7 convalescent emc_prnt' = 'B.1.1.7',
  'B.1.1.7 convalescent kcl' = 'B.1.1.7',
  'B.1.1.7 convalescent innsbruck' ='B.1.1.7',
  'B.1.1.7 convalescent fda' = 'B.1.1.7',
  'B.1.1.7 convalescent geneva' = 'B.1.1.7',
  'B.1.1.7 convalescent amc' = 'B.1.1.7',
  'B.1.1.7 convalescent charite' = 'B.1.1.7',
  'B.1.1.7 convalescent madison_frnt' = 'B.1.1.7',
  'B.1.1.7 convalescent emc_calu' = 'B.1.1.7',
  'B.1.1.7 convalescent emc_vero' = 'B.1.1.7',
  'BA.1 convalescent duke' = 'BA.1',
  'BA.1 convalescent maryland' = 'BA.1',
  'BA.1 convalescent galveston' = 'BA.1',
  'BA.1 convalescent emory' = 'BA.1',
  'BA.1 convalescent madison_pooled' = 'BA.1',
  'BA.1 convalescent madison_unpooled' = 'BA.1',
  'BA.1 convalescent st_louis' = 'BA.1',
  'BA.1 convalescent oxford' = 'BA.1',
  'BA.1 convalescent mt_sinai_mouse' = 'BA.1',
  'BA.1 convalescent mt_sinai_human' = 'BA.1',
  'BA.1 convalescent emc_prnt' = 'BA.1',
  'BA.1 convalescent kcl' = 'BA.1',
  'BA.1 convalescent innsbruck' ='BA.1',
  'BA.1 convalescent fda' = 'BA.1',
  'BA.1 convalescent geneva' = 'BA.1',
  'BA.1 convalescent amc' = 'BA.1',
  'BA.1 convalescent charite' = 'BA.1',
  'BA.1 convalescent madison_frnt' = 'BA.1',
  'BA.1 convalescent emc_calu' = 'BA.1',
  'BA.1 convalescent emc_vero' = 'BA.1'
)

# Read in data
maps <- load_maps_for_comparison('', sr_group_standardisation = 'data/metadata/sr_group_standardisation_WT_conv_as_D614G_as614D.csv')

slope_sr_groups <- c('D614G convalescent', 'mRNA-1273', 'B.1.1.7 convalescent', 'B.1.351 convalescent', 
                     'P.1 convalescent', 'B.1.617.2 convalescent', 'BA.1 convalescent')

# Arrange the data
data_orig <- arrange_data(maps, dilution_stepsize = 0)
data <- subsetMap(data_orig, sera = srGroups(data_orig) %in% slope_sr_groups, antigens = agNames(data_orig) %in% duplicated_ags)

srGroups(data) <- factor(srGroups(data), levels = slope_sr_groups)

# Run the model
options(mc.cores = parallel::detectCores())

# Build up titer information
titer_layers <- Racmacs:::titerTableLayers(data)

upper_logtiter_lims <- Racmacs:::logtiterTableLayers(data)
lower_logtiter_lims <- Racmacs:::logtiterTableLayers(data)

for (map_ in layerNames(data)) {
  for (serum_group in slope_sr_groups) {
    # get homologous ag
    homologous_ag <- unname(homologous_ags[paste(serum_group, map_)])
    
    for (ag in agNames(data)) {
      diff_lims <- titertools:::calc_titer_diff_lims(
        titers1 = titer_layers[[map_]][match(homologous_ag, agNames(data)), which(srGroups(data) %in% c(serum_group))],
        titers2 = titer_layers[[map_]][match(ag, agNames(data)), which(srGroups(data) %in% c(serum_group))],
        dilution_stepsize = dilutionStepsize(maps[[map_]])
      )
      upper_logtiter_lims[[map_]][match(ag, agNames(data)), which(srGroups(data) %in% c(serum_group))] <- diff_lims$max_diffs
      lower_logtiter_lims[[map_]][match(ag, agNames(data)), which(srGroups(data) %in% c(serum_group))] <- diff_lims$min_diffs
    }
  }
}

# Replace NAs with 0s
upper_logtiter_lims <- lapply(upper_logtiter_lims, \(upper_lims) {
  upper_lims[is.na(upper_lims)] <- 0
  upper_lims
})

lower_logtiter_lims <- lapply(lower_logtiter_lims, \(lower_lims) {
  lower_lims[is.na(lower_lims)] <- 0
  lower_lims
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
  titertypes = simplify2array(Racmacs:::titertypesTableLayers(data)),
  lower_logtiter_lims = simplify2array(lower_logtiter_lims),
  upper_logtiter_lims = simplify2array(upper_logtiter_lims)
)

# Decide on starting conditions for the model
data_init <- list(
  logtiter_error_sigma = rep(0.8, numLayers(data))
)

# Fetch the model
mod <- cmdstan_model(stan_file = 'code/stan_models/slope_calculation_by_animal_assay.stan')

# Sample
mod_optimized <- mod$sample(
  data = data_list,
  init = list(data_init, data_init, data_init, data_init),
  chains = 4, 
  refresh = 500,
  iter_sampling  = 10000,
  iter_warmup = 1000
)

# Draw from the posterior
mod_optimized_draws_animal <- mod_optimized$draws(format = 'df', variables = c('animal_slope_effect'))
mod_optimized_draws_assay <- mod_optimized$draws(format = 'df', variables = c('assay_slope_effect'))

# Re-name columns
colnames(mod_optimized_draws_animal) <- c('human', 'hamster', 'mouse', '.chain', '.iteration', '.draw')
colnames(mod_optimized_draws_assay) <- c('FRNT', 'PV-neut', 'VSV-PV-neut', 'PRNT', 'Microneut', 'CPE', '.chain', '.iteration', '.draw')

# Save the draws
saveRDS(mod_optimized_draws_animal, 'data/titer_analyses/foldchange/slope_calculation_by_animal_assay/slope_calculation_draws_animal.rds')
saveRDS(mod_optimized_draws_assay, 'data/titer_analyses/foldchange/slope_calculation_by_animal_assay/slope_calculation_draws_assay.rds')

# Check convergence
summary <- mod_optimized$summary()
summary
# # A tibble: 189 × 10
# variable                          mean   median      sd     mad          q5         q95  rhat ess_b…¹ ess_t…²
# <chr>                            <dbl>    <dbl>   <dbl>   <dbl>       <dbl>       <dbl> <dbl>   <dbl>   <dbl>
#   1 lp__                     -10855.       -1.09e+4 12.8    12.8    -10877.        -1.08e+4  1.00   1938.   6294.
# 2 animal_slope_effect[1]     0.179     0.144  0.129  0.0603     0.0805     0.388  1.00    2724.    2489.
# 3 animal_slope_effect[2]     0.141     0.114  0.102  0.0477     0.0633     0.305  1.00    2727.    2500.
# 4 animal_slope_effect[3]     0.139     0.112  0.101  0.0476     0.0617     0.302  1.00    2747.    2587.
# 5 assay_slope_effect[1]      6.82      6.62   2.85   2.91       2.46      11.8    1.00    2795.    2489.
# 6 assay_slope_effect[2]      9.06      8.80   3.79   3.86       3.26      15.7    1.00    2799.    2457.
# 7 assay_slope_effect[3]      6.50      6.32   2.72   2.77       2.34      11.3    1.00    2804.    2469.
# 8 assay_slope_effect[4]      6.92      6.72   2.90   2.95       2.50      12.0    1.00    2799.    2490.
# 9 assay_slope_effect[5]      7.65      7.42   3.21   3.27       2.75      13.3    1.00    2817.    2484.
# 10 assay_slope_effect[6]      3.79      3.65   1.62   1.64       1.35       6.65   1.00    2850.    2523.

# Look at traces
gp <- mcmc_trace(mod_optimized$draws(variables = c('animal_slope_effect', 'assay_slope_effect')))

ggsave('data/titer_analyses/foldchange/slope_calculation_by_animal_assay/animal_assay_slope_effects_traceplots.png', gp, width = 15, height = 10)

