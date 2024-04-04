
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
  titertypes = simplify2array(Racmacs:::titertypesTableLayers(data)),
  lower_logtiter_lims = simplify2array(lower_logtiter_lims),
  upper_logtiter_lims = simplify2array(upper_logtiter_lims)
)

# Decide on starting conditions for the model
data_init <- list(
  logtiter_error_sigma = rep(0.8, numLayers(data)),
  dataset_slope_effect = rep(0, numLayers(data))
  # sr_group_slopes = rep(0.9, numSeraGroups(data)),
  # ag_folddrops = rep(1, numAntigens(data))
)

# Fetch the model
mod <- cmdstan_model(stan_file = 'code/stan_models/slope_calculation.stan')

# Sample
mod_optimized <- mod$sample(
  data = data_list,
  init = list(data_init, data_init, data_init, data_init),
  chains = 4, 
  refresh = 500,
  iter_sampling  = 10000,
  iter_warmup = 1000
)

# Draw from the posterior, map slope effect
mod_optimized_draws <- mod_optimized$draws(format = 'df', variables = c('dataset_slope_effect'))

# Re-name columns
colnames(mod_optimized_draws) <- c('duke', 'maryland', 'galveston', 'emory',
                                   'madison_pooled', 'madison_unpooled', 'st_louis',
                                   'oxford', 'mt_sinai_human', 'emc_prnt',
                                   'innsbruck', 'charite', 'madison_frnt', 'fda', 'geneva',
                                   'amc', 'emc_calu', 'emc_vero', '.chain', '.iteration', '.draw')

# Save the draws
saveRDS(mod_optimized_draws, 'data/titer_analyses/foldchange/slope_calculation/slope_calculation_dataset_slope_effect_draws.rds')

# Draw from the posterior, ag_folddrops
mod_optimized_draws_foldchange <- mod_optimized$draws(format = 'df', variables = c('ag_folddrops'))

# Save the draws
saveRDS(mod_optimized_draws_foldchange, 'data/titer_analyses/foldchange/slope_calculation/slope_calculation_ag_folddrops_draws.rds')

# Check convergence
summary <- mod_optimized$summary()
print(summary, n=200)
# # A tibble: 198 × 10
# variable                         mean       median      sd     mad        q5      q95  rhat ess_bulk ess_tail
# <chr>                           <dbl>        <dbl>   <dbl>   <dbl>     <dbl>    <dbl> <dbl>    <dbl>    <dbl>
#   1 lp__                     -9332.       -9332.       13.4    13.4    -9355.    -9.31e+3  1.00    1534.    5391.
# 2 dataset_slope_effect[1]      1.45         1.45      0.0904  0.0898     1.31   1.61e+0  1.00     808.    1656.
# 3 dataset_slope_effect[2]      0.397        0.395     0.0947  0.0929     0.244  5.55e-1  1.00    9948.   16985.
# 4 dataset_slope_effect[3]      0.848        0.845     0.0750  0.0742     0.730  9.76e-1  1.00    1535.    4097.
# 5 dataset_slope_effect[4]      1.04         1.03      0.0656  0.0656     0.935  1.15e+0  1.00     818.    1639.
# 6 dataset_slope_effect[5]      0.502        0.501     0.0956  0.0946     0.348  6.63e-1  1.00    7108.   17266.
# 7 dataset_slope_effect[6]      0.449        0.447     0.0563  0.0557     0.359  5.45e-1  1.00    3060.    7425.
# 8 dataset_slope_effect[7]      0.840        0.837     0.0817  0.0803     0.710  9.80e-1  1.00    1791.    4585.
# 9 dataset_slope_effect[8]      1.21         1.21      0.0820  0.0813     1.09   1.36e+0  1.00     909.    1874.
# 10 dataset_slope_effect[9]      1.20         1.19      0.0901  0.0893     1.06   1.35e+0  1.00    1101.    2783.
# 11 dataset_slope_effect[10]     0.778        0.775     0.0588  0.0585     0.686  8.80e-1  1.00    1197.    2835.
# 12 dataset_slope_effect[11]     0.842        0.839     0.0582  0.0577     0.751  9.43e-1  1.00     996.    2208.
# 13 dataset_slope_effect[12]     0.786        0.783     0.0582  0.0582     0.696  8.86e-1  1.00    1121.    2783.
# 14 dataset_slope_effect[13]     0.789        0.786     0.0590  0.0586     0.698  8.91e-1  1.00    1182.    2921.
# 15 dataset_slope_effect[14]     1.34         1.33      0.0838  0.0830     1.21   1.49e+0  1.00     830.    1606.
# 16 dataset_slope_effect[15]     1.04         1.04      0.0687  0.0687     0.932  1.16e+0  1.00     905.    1880.
# 17 dataset_slope_effect[16]     0.978        0.975     0.0658  0.0649     0.876  1.09e+0  1.00     935.    2065.
# 18 dataset_slope_effect[17]     0.763        0.761     0.0539  0.0539     0.679  8.56e-1  1.00    1044.    2306.
# 19 dataset_slope_effect[18]     0.722        0.719     0.0566  0.0563     0.633  8.20e-1  1.00    1249.    3034.
# 20 ag_folddrops[1,1]           -0.287       -0.286     0.0967  0.0952    -0.447 -1.29e-1  1.00   22476.   24144.
# 21 ag_folddrops[2,1]           -0.848       -0.845     0.108   0.107     -1.03  -6.78e-1  1.00    3309.    9180.
# 22 ag_folddrops[3,1]           -2.37        -2.36      0.192   0.194     -2.69  -2.06e+0  1.00    1321.    3334.
# 23 ag_folddrops[4,1]           -1.62        -1.62      0.139   0.139     -1.86  -1.40e+0  1.00    1506.    3976.
# 24 ag_folddrops[5,1]           -2.78        -2.78      0.197   0.196     -3.11  -2.46e+0  1.00    1015.    2167.
# 25 ag_folddrops[6,1]           -0.698       -0.695     0.152   0.151     -0.953 -4.51e-1  1.00    8558.   16671.
# 26 ag_folddrops[7,1]           -1.61        -1.61      0.178   0.178     -1.91  -1.33e+0  1.00    2562.    7905.
# 27 ag_folddrops[8,1]           -2.25        -2.24      0.223   0.222     -2.62  -1.89e+0  1.00    1982.    6684.
# 28 ag_folddrops[9,1]           -1.17        -1.17      0.116   0.115     -1.37  -9.85e-1  1.00    1924.    5836.
# 29 ag_folddrops[10,1]          -1.00        -0.996     0.187   0.184     -1.31  -7.01e-1  1.00    6521.   14898.
# 30 ag_folddrops[11,1]          -2.64        -2.64      0.227   0.229     -3.03  -2.28e+0  1.00    1473.    4713.
# 31 ag_folddrops[12,1]          -1.02        -1.01      0.160   0.158     -1.29  -7.57e-1  1.00    5160.   13212.
# 32 ag_folddrops[13,1]          -5.94        -5.93      0.401   0.403     -6.61  -5.29e+0  1.00     954.    2101.
# 33 ag_folddrops[14,1]          -6.97        -6.94      0.644   0.638     -8.06  -5.94e+0  1.00    1975.    5864.
# 34 ag_folddrops[15,1]          -4.33        -4.32      0.325   0.325     -4.88  -3.81e+0  1.00    1180.    2654.
# 35 ag_folddrops[16,1]          -5.53        -5.51      0.593   0.592     -6.55  -4.59e+0  1.00    2617.    8082.
# 36 ag_folddrops[17,1]          -5.17        -5.15      0.484   0.479     -6.00  -4.40e+0  1.00    1896.    5900.
# 37 ag_folddrops[18,1]          -0.140       -0.139     0.126   0.124     -0.348  6.65e-2  1.00   48433.   26352.
# 38 ag_folddrops[19,1]          -1.55        -1.55      0.218   0.218     -1.92  -1.20e+0  1.00    3812.   10509.
# 39 ag_folddrops[20,1]           0.0175       0.0179    0.203   0.201     -0.317  3.53e-1  1.00   53173.   27621.
# 40 ag_folddrops[21,1]           0.143        0.143     0.204   0.203     -0.191  4.79e-1  1.00   46067.   27940.
# 41 ag_folddrops[22,1]          -2.20        -2.19      0.192   0.192     -2.52  -1.89e+0  1.00    1500.    4624.
# 42 ag_folddrops[23,1]           0.695        0.693     1.24    1.23      -1.33   2.73e+0  1.00   46207.   29197.
# 43 ag_folddrops[1,2]           -0.000843    -0.000892  0.122   0.121     -0.200  2.00e-1  1.00   49725.   30573.
# 44 ag_folddrops[2,2]           -0.701       -0.700     0.123   0.123     -0.907 -5.02e-1  1.00    5771.   19785.
# 45 ag_folddrops[3,2]           -1.63        -1.62      0.140   0.140     -1.86  -1.40e+0  1.00    1415.    3660.
# 46 ag_folddrops[4,2]           -1.39        -1.38      0.143   0.144     -1.63  -1.16e+0  1.00    2131.    6932.
# 47 ag_folddrops[5,2]           -2.11        -2.11      0.166   0.166     -2.39  -1.84e+0  1.00    1217.    3530.
# 48 ag_folddrops[6,2]           -1.18        -1.17      0.173   0.172     -1.47  -9.00e-1  1.00    4216.   11789.
# 49 ag_folddrops[7,2]           -0.877       -0.873     0.118   0.118     -1.08  -6.90e-1  1.00    3412.   10550.
# 50 ag_folddrops[8,2]           -2.08        -2.08      0.200   0.199     -2.42  -1.76e+0  1.00    1810.    6532.
# 51 ag_folddrops[9,2]           -1.28        -1.27      0.117   0.117     -1.47  -1.09e+0  1.00    1621.    5193.
# 52 ag_folddrops[10,2]          -0.974       -0.970     0.188   0.188     -1.29  -6.72e-1  1.00    7756.   18136.
# 53 ag_folddrops[11,2]          -2.06        -2.05      0.220   0.220     -2.43  -1.71e+0  1.00    2192.    6596.
# 54 ag_folddrops[12,2]          -1.44        -1.44      0.140   0.140     -1.68  -1.22e+0  1.00    1847.    5959.
# 55 ag_folddrops[13,2]          -5.55        -5.54      0.388   0.389     -6.20  -4.92e+0  1.00    1029.    2425.
# 56 ag_folddrops[14,2]          -3.59        -3.59      0.622   0.619     -4.63  -2.58e+0  1.00    6625.   17047.
# 57 ag_folddrops[15,2]          -3.87        -3.85      0.644   0.639     -4.96  -2.83e+0  1.00    6176.   18736.
# 58 ag_folddrops[16,2]          -3.88        -3.87      0.627   0.622     -4.93  -2.86e+0  1.00    5822.   15102.
# 59 ag_folddrops[17,2]          -4.58        -4.57      0.505   0.508     -5.43  -3.78e+0  1.00    2501.    8310.
# 60 ag_folddrops[18,2]          -0.0523      -0.0513    0.0930  0.0925    -0.206  9.99e-2  1.00   53027.   26575.
# 61 ag_folddrops[19,2]          -0.623       -0.620     0.202   0.202     -0.958 -2.97e-1  1.00   20842.   24437.
# 62 ag_folddrops[20,2]          -0.161       -0.160     0.106   0.105     -0.335  1.16e-2  1.00   42885.   26254.
# 63 ag_folddrops[21,2]          -0.781       -0.778     0.118   0.118     -0.979 -5.92e-1  1.00    4570.   14691.
# 64 ag_folddrops[22,2]          -1.93        -1.92      0.178   0.179     -2.23  -1.64e+0  1.00    1720.    5813.
# 65 ag_folddrops[23,2]          -1.00        -0.993     3.02    3.03      -5.98   3.98e+0  1.00   48516.   30029.
# 66 ag_folddrops[1,3]           -1.10        -1.10      0.161   0.161     -1.37  -8.40e-1  1.00    4191.   11238.
# 67 ag_folddrops[2,3]           -0.00230     -0.00182   0.147   0.146     -0.245  2.40e-1  1.00   51941.   28891.
# 68 ag_folddrops[3,3]           -0.972       -0.969     0.198   0.197     -1.30  -6.54e-1  1.00    8085.   16873.
# 69 ag_folddrops[4,3]           -1.55        -1.54      0.184   0.185     -1.86  -1.25e+0  1.00    2762.    8319.
# 70 ag_folddrops[5,3]           -2.13        -2.13      0.201   0.201     -2.47  -1.81e+0  1.00    1814.    6217.
# 71 ag_folddrops[6,3]           -1.63        -1.62      0.298   0.296     -2.13  -1.14e+0  1.00    7032.   18053.
# 72 ag_folddrops[7,3]           -3.69        -3.68      0.387   0.387     -4.35  -3.08e+0  1.00    2390.    8170.
# 73 ag_folddrops[8,3]           -2.97        -2.96      0.357   0.357     -3.58  -2.40e+0  1.00    3114.   10352.
# 74 ag_folddrops[9,3]           -1.68        -1.67      0.183   0.182     -1.99  -1.39e+0  1.00    2386.    7535.
# 75 ag_folddrops[10,3]          -1.92        -1.91      0.352   0.350     -2.51  -1.35e+0  1.00    6899.   15660.
# 76 ag_folddrops[11,3]          -2.12        -2.11      0.299   0.299     -2.62  -1.64e+0  1.00    4083.   13095.
# 77 ag_folddrops[12,3]          -1.47        -1.47      0.287   0.287     -1.95  -1.01e+0  1.00    7808.   18131.
# 78 ag_folddrops[13,3]          -6.67        -6.66      0.486   0.484     -7.49  -5.88e+0  1.00    1128.    3009.
# 79 ag_folddrops[14,3]          -7.00        -6.82      1.36    1.26      -9.52  -5.12e+0  1.00   15086.   20677.
# 80 ag_folddrops[15,3]          -6.00        -5.99      0.470   0.471     -6.80  -5.25e+0  1.00    1280.    3554.
# 81 ag_folddrops[16,3]          -6.99        -6.81      1.35    1.25      -9.52  -5.10e+0  1.00   14252.   20265.
# 82 ag_folddrops[17,3]          -5.22        -5.20      0.546   0.540     -6.15  -4.36e+0  1.00    2410.    8063.
# 83 ag_folddrops[18,3]           0.144        0.143     0.229   0.229     -0.230  5.19e-1  1.00   48037.   27998.
# 84 ag_folddrops[19,3]          -3.67        -3.66      0.675   0.673     -4.80  -2.58e+0  1.00    8013.   18958.
# 85 ag_folddrops[20,3]          -0.992       -0.992     3.00    2.99      -5.91   3.91e+0  1.00   57034.   27841.
# 86 ag_folddrops[21,3]          -3.02        -3.02      2.30    2.29      -6.82   8.04e-1  1.00   45831.   28374.
# 87 ag_folddrops[22,3]          -2.61        -2.60      0.278   0.275     -3.08  -2.17e+0  1.00    2311.    7388.
# 88 ag_folddrops[23,3]          -3.03        -3.03      2.31    2.30      -6.77   7.62e-1  1.00   47681.   28100.
# 89 ag_folddrops[1,4]           -2.68        -2.67      0.205   0.201     -3.02  -2.35e+0  1.00    1198.    2998.
# 90 ag_folddrops[2,4]           -1.24        -1.24      0.150   0.149     -1.50  -1.00e+0  1.00    2810.    8195.
# 91 ag_folddrops[3,4]           -0.916       -0.913     0.136   0.135     -1.14  -6.97e-1  1.00    4098.   10775.
# 92 ag_folddrops[4,4]           -0.231       -0.229     0.132   0.132     -0.448 -1.64e-2  1.00   33361.   23908.
# 93 ag_folddrops[5,4]           -0.00161     -0.00128   0.103   0.102     -0.171  1.67e-1  1.00   51249.   27936.
# 94 ag_folddrops[6,4]           -2.12        -2.12      0.216   0.215     -2.49  -1.78e+0  1.00    2148.    6869.
# 95 ag_folddrops[7,4]           -1.79        -1.79      0.187   0.186     -2.11  -1.49e+0  1.00    2227.    7075.
# 96 ag_folddrops[8,4]           -2.64        -2.64      0.253   0.252     -3.07  -2.24e+0  1.00    1895.    5795.
# 97 ag_folddrops[9,4]           -3.09        -3.08      0.228   0.228     -3.47  -2.72e+0  1.00    1099.    2954.
# 98 ag_folddrops[10,4]          -2.14        -2.14      0.243   0.243     -2.55  -1.75e+0  1.00    2541.    7825.
# 99 ag_folddrops[11,4]          -0.794       -0.791     0.169   0.168     -1.08  -5.21e-1  1.00    9387.   17807.
# 100 ag_folddrops[12,4]          -2.63        -2.63      0.225   0.227     -3.01  -2.27e+0  1.00    1499.    4368.
# 101 ag_folddrops[13,4]          -4.43        -4.42      0.333   0.332     -4.99  -3.89e+0  1.00    1173.    3028.
# 102 ag_folddrops[14,4]          -3.07        -3.06      0.385   0.382     -3.72  -2.45e+0  1.00    3332.   10194.
# 103 ag_folddrops[15,4]          -3.76        -3.75      0.361   0.362     -4.36  -3.18e+0  1.00    1931.    6389.
# 104 ag_folddrops[16,4]          -3.65        -3.63      0.475   0.471     -4.44  -2.89e+0  1.00    3647.   11494.
# 105 ag_folddrops[17,4]          -4.14        -4.13      0.389   0.386     -4.79  -3.52e+0  1.00    1941.    5987.
# 106 ag_folddrops[18,4]          -2.67        -2.67      0.227   0.227     -3.06  -2.31e+0  1.00    1440.    4317.
# 107 ag_folddrops[19,4]          -0.826       -0.823     0.219   0.217     -1.19  -4.73e-1  1.00   11350.   20810.
# 108 ag_folddrops[20,4]          -2.94        -2.93      0.276   0.276     -3.40  -2.50e+0  1.00    1733.    5255.
# 109 ag_folddrops[21,4]          -2.77        -2.76      0.260   0.259     -3.21  -2.36e+0  1.00    1720.    5525.
# 110 ag_folddrops[22,4]          -1.78        -1.78      0.202   0.201     -2.13  -1.46e+0  1.00    2643.    8901.
# 111 ag_folddrops[23,4]          -3.11        -3.09      1.34    1.33      -5.36  -9.44e-1  1.00   36903.   26486.
# 112 ag_folddrops[1,5]           -2.65        -2.64      0.260   0.260     -3.09  -2.24e+0  1.00    2010.    7106.
# 113 ag_folddrops[2,5]           -2.06        -2.05      0.231   0.231     -2.45  -1.69e+0  1.00    2660.    7206.
# 114 ag_folddrops[3,5]           -2.13        -2.12      0.252   0.250     -2.56  -1.72e+0  1.00    2907.   10313.
# 115 ag_folddrops[4,5]           -0.00302     -0.00335   0.179   0.179     -0.295  2.94e-1  1.00   54656.   28229.
# 116 ag_folddrops[5,5]           -0.797       -0.795     0.188   0.187     -1.11  -4.94e-1  1.00   10728.   18056.
# 117 ag_folddrops[6,5]           -3.63        -3.62      0.372   0.371     -4.27  -3.05e+0  1.00    2244.    7382.
# 118 ag_folddrops[7,5]           -1.74        -1.74      0.317   0.317     -2.27  -1.23e+0  1.00    6196.   14442.
# 119 ag_folddrops[8,5]           -3.24        -3.23      0.406   0.407     -3.93  -2.60e+0  1.00    3319.    9376.
# 120 ag_folddrops[9,5]           -4.18        -4.18      0.338   0.340     -4.75  -3.64e+0  1.00    1363.    3383.
# 121 ag_folddrops[10,5]          -2.78        -2.77      0.385   0.387     -3.43  -2.17e+0  1.00    3999.   11714.
# 122 ag_folddrops[11,5]          -1.25        -1.24      0.323   0.321     -1.79  -7.26e-1  1.00   13401.   21939.
# 123 ag_folddrops[12,5]          -3.12        -3.11      0.379   0.377     -3.76  -2.51e+0  1.00    3149.   10244.
# 124 ag_folddrops[13,5]          -4.67        -4.66      0.418   0.418     -5.37  -4.00e+0  1.00    1721.    4962.
# 125 ag_folddrops[14,5]          -3.65        -3.64      0.715   0.709     -4.85  -2.50e+0  1.00    8684.   19149.
# 126 ag_folddrops[15,5]          -4.20        -4.19      0.440   0.437     -4.95  -3.50e+0  1.00    2430.    8448.
# 127 ag_folddrops[16,5]          -3.49        -3.48      0.703   0.694     -4.68  -2.37e+0  1.00    9678.   20119.
# 128 ag_folddrops[17,5]          -3.19        -3.18      0.689   0.687     -4.34  -2.08e+0  1.00   11168.   19394.
# 129 ag_folddrops[18,5]          -2.70        -2.69      0.326   0.328     -3.25  -2.18e+0  1.00    2838.    9882.
# 130 ag_folddrops[19,5]          -1.37        -1.37      0.738   0.727     -2.58  -1.63e-1  1.00   31532.   23822.
# 131 ag_folddrops[20,5]          -0.992       -0.994     3.04    3.01      -6.00   4.02e+0  1.00   57210.   28729.
# 132 ag_folddrops[21,5]          -3.03        -3.04      2.32    2.28      -6.84   7.92e-1  1.00   46848.   27528.
# 133 ag_folddrops[22,5]          -2.73        -2.72      0.319   0.317     -3.27  -2.22e+0  1.00    2767.    9671.
# 134 ag_folddrops[23,5]          -3.02        -3.03      2.30    2.28      -6.80   7.90e-1  1.00   49839.   29524.
# 135 ag_folddrops[1,6]           -2.78        -2.78      0.225   0.223     -3.16  -2.42e+0  1.00    1357.    3562.
# 136 ag_folddrops[2,6]           -2.93        -2.93      0.236   0.236     -3.33  -2.55e+0  1.00    1357.    3962.
# 137 ag_folddrops[3,6]           -4.73        -4.72      0.419   0.419     -5.44  -4.06e+0  1.00    1656.    4401.
# 138 ag_folddrops[4,6]           -4.17        -4.16      0.306   0.307     -4.68  -3.68e+0  1.00    1124.    2995.
# 139 ag_folddrops[5,6]           -4.56        -4.55      0.318   0.318     -5.09  -4.04e+0  1.00    1000.    2390.
# 140 ag_folddrops[6,6]           -2.52        -2.52      0.246   0.245     -2.94  -2.13e+0  1.00    1964.    6398.
# 141 ag_folddrops[7,6]           -4.33        -4.33      0.331   0.332     -4.89  -3.80e+0  1.00    1215.    3467.
# 142 ag_folddrops[8,6]           -3.89        -3.88      0.326   0.328     -4.44  -3.37e+0  1.00    1411.    4577.
# 143 ag_folddrops[9,6]           -0.00125     -0.00129   0.145   0.144     -0.237  2.37e-1  1.00   50296.   28173.
# 144 ag_folddrops[10,6]          -0.156       -0.157     0.310   0.309     -0.665  3.53e-1  1.00   55071.   29024.
# 145 ag_folddrops[11,6]          -4.80        -4.79      0.351   0.350     -5.38  -4.23e+0  1.00    1125.    2788.
# 146 ag_folddrops[12,6]          -3.13        -3.12      0.269   0.272     -3.58  -2.70e+0  1.00    1554.    4979.
# 147 ag_folddrops[13,6]          -5.09        -5.09      0.342   0.345     -5.66  -4.53e+0  1.00     963.    2203.
# 148 ag_folddrops[14,6]          -4.70        -4.69      0.378   0.378     -5.33  -4.09e+0  1.00    1372.    3901.
# 149 ag_folddrops[15,6]          -4.62        -4.62      0.351   0.352     -5.21  -4.06e+0  1.00    1226.    3552.
# 150 ag_folddrops[16,6]          -5.10        -5.09      0.506   0.500     -5.95  -4.28e+0  1.00    2133.    7068.
# 151 ag_folddrops[17,6]          -4.72        -4.71      0.393   0.393     -5.38  -4.09e+0  1.00    1470.    4179.
# 152 ag_folddrops[18,6]          -3.47        -3.46      0.362   0.359     -4.08  -2.89e+0  1.00    2286.    7525.
# 153 ag_folddrops[19,6]          -3.54        -3.54      0.417   0.414     -4.25  -2.88e+0  1.00    2850.   10145.
# 154 ag_folddrops[20,6]          -3.30        -3.29      0.419   0.421     -4.00  -2.63e+0  1.00    3347.   11411.
# 155 ag_folddrops[21,6]          -2.95        -2.94      0.400   0.397     -3.63  -2.32e+0  1.00    3842.    9303.
# 156 ag_folddrops[22,6]          -3.57        -3.57      0.356   0.355     -4.17  -3.00e+0  1.00    2015.    6379.
# 157 ag_folddrops[23,6]           0.157        0.152     1.33    1.33      -2.02   2.34e+0  1.00   51826.   28545.
# 158 ag_folddrops[1,7]           -5.40        -5.39      0.492   0.488     -6.22  -4.61e+0  1.00    1778.    5312.
# 159 ag_folddrops[2,7]           -5.08        -5.07      0.491   0.493     -5.91  -4.30e+0  1.00    2012.    6701.
# 160 ag_folddrops[3,7]           -4.18        -4.14      0.874   0.872     -5.66  -2.80e+0  1.00    9795.   20988.
# 161 ag_folddrops[4,7]           -3.92        -3.91      0.562   0.566     -4.87  -3.02e+0  1.00    4757.   13455.
# 162 ag_folddrops[5,7]           -5.01        -5.00      0.457   0.454     -5.79  -4.28e+0  1.00    1781.    5527.
# 163 ag_folddrops[6,7]           -3.95        -3.91      0.937   0.936     -5.56  -2.48e+0  1.00   16192.   23861.
# 164 ag_folddrops[7,7]           -3.85        -3.81      0.930   0.920     -5.45  -2.39e+0  1.00   16915.   24713.
# 165 ag_folddrops[8,7]           -2.84        -2.83      0.554   0.552     -3.78  -1.95e+0  1.00    7618.   20955.
# 166 ag_folddrops[9,7]           -4.82        -4.81      0.420   0.424     -5.53  -4.15e+0  1.00    1655.    5352.
# 167 ag_folddrops[10,7]          -0.994       -0.984     2.99    3.01      -5.90   3.87e+0  1.00   56993.   29665.
# 168 ag_folddrops[11,7]          -4.36        -4.34      0.635   0.634     -5.43  -3.34e+0  1.00    4996.   11836.
# 169 ag_folddrops[12,7]          -2.86        -2.84      0.664   0.664     -3.97  -1.80e+0  1.00   13611.   18663.
# 170 ag_folddrops[13,7]          -0.00835     -0.00752   0.247   0.242     -0.412  3.95e-1  1.00   51578.   27862.
# 171 ag_folddrops[14,7]          -0.269       -0.269     0.342   0.341     -0.835  2.91e-1  1.00   47522.   29538.
# 172 ag_folddrops[15,7]          -1.94        -1.93      0.303   0.302     -2.45  -1.45e+0  1.00    5052.   13342.
# 173 ag_folddrops[16,7]          -0.693       -0.689     0.548   0.545     -1.60   2.02e-1  1.00   45965.   27178.
# 174 ag_folddrops[17,7]          -3.54        -3.53      0.467   0.468     -4.33  -2.80e+0  1.00    3596.   12313.
# 175 ag_folddrops[18,7]          -6.73        -6.71      0.603   0.602     -7.76  -5.77e+0  1.00    1783.    5476.
# 176 ag_folddrops[19,7]          -0.977       -0.984     2.99    2.97      -5.92   3.97e+0  1.00   51537.   29871.
# 177 ag_folddrops[20,7]          -0.965       -0.976     3.01    3.01      -5.91   4.00e+0  1.00   48810.   28386.
# 178 ag_folddrops[21,7]          -0.979       -0.987     3.02    3.01      -5.96   3.99e+0  1.00   49326.   29518.
# 179 ag_folddrops[22,7]          -3.80        -3.78      1.00    0.988     -5.49  -2.18e+0  1.00   19012.   23637.
# 180 ag_folddrops[23,7]          -1.00        -0.998     3.00    2.99      -5.93   3.92e+0  1.00   52807.   29500.
# 181 logtiter_error_sigma[1]      1.72         1.72      0.0414  0.0413     1.65   1.79e+0  1.00   41024.   30554.
# 182 logtiter_error_sigma[2]      1.53         1.51      0.201   0.195      1.23   1.89e+0  1.00   44871.   27313.
# 183 logtiter_error_sigma[3]      0.929        0.923     0.0835  0.0824     0.802  1.08e+0  1.00   40432.   29782.
# 184 logtiter_error_sigma[4]      0.850        0.849     0.0246  0.0246     0.811  8.92e-1  1.00   39639.   30978.
# 185 logtiter_error_sigma[5]      1.70         1.69      0.175   0.173      1.43   2.00e+0  1.00   41470.   28511.
# 186 logtiter_error_sigma[6]      1.10         1.10      0.0921  0.0912     0.958  1.26e+0  1.00   50773.   29049.
# 187 logtiter_error_sigma[7]      1.17         1.16      0.0751  0.0740     1.05   1.30e+0  1.00   44294.   28713.
# 188 logtiter_error_sigma[8]      1.17         1.17      0.0377  0.0375     1.11   1.23e+0  1.00   38229.   30643.
# 189 logtiter_error_sigma[9]      0.767        0.765     0.0405  0.0404     0.703  8.36e-1  1.00   37129.   31083.
# 190 logtiter_error_sigma[10]     1.58         1.58      0.0843  0.0834     1.45   1.73e+0  1.00   44570.   28392.
# 191 logtiter_error_sigma[11]     1.54         1.53      0.0660  0.0657     1.43   1.65e+0  1.00   45619.   28472.
# 192 logtiter_error_sigma[12]     1.18         1.17      0.0949  0.0935     1.03   1.34e+0  1.00   45005.   28858.
# 193 logtiter_error_sigma[13]     1.13         1.12      0.101   0.0988     0.973  1.30e+0  1.00   45568.   28407.
# 194 logtiter_error_sigma[14]     1.65         1.65      0.0567  0.0567     1.56   1.74e+0  1.00   44735.   30795.
# 195 logtiter_error_sigma[15]     1.60         1.60      0.0555  0.0555     1.51   1.69e+0  1.00   42190.   31677.
# 196 logtiter_error_sigma[16]     1.37         1.37      0.0659  0.0653     1.27   1.49e+0  1.00   40805.   28067.
# 197 logtiter_error_sigma[17]     1.05         1.05      0.0687  0.0685     0.942  1.17e+0  1.00   43241.   29671.
# 198 logtiter_error_sigma[18]     1.42         1.41      0.0886  0.0879     1.28   1.57e+0  1.00   44143.   30117.


# Look at traces
gp <- mcmc_trace(mod_optimized$draws(variables = c('dataset_slope_effect')))

ggsave('data/titer_analyses/foldchange/slope_calculation/slope_calculation_dataset_slope_effect_traceplots.png', gp, width = 15, height = 10)
