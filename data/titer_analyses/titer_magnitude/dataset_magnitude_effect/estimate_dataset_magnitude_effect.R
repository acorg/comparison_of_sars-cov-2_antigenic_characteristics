
rm(list = ls())

library(Racmacs)
library(tidyverse)
library(cmdstanr)
library(titertools)

source('code/data_generation/load_maps_for_comparison.R')
source('code/functions/gmts.R')
source('code/functions/map_longinfo.R')

# Read in data
maps <- load_maps_for_comparison('', sr_group_standardisation = 'data/metadata/sr_group_standardisation_WT_conv_as_D614G_as614D.csv')

# Arrange the data
merged_map <- arrange_data(maps, dilution_stepsize = 0)

# Run the model
options(mc.cores = parallel::detectCores())

# Replace logtiter table layer NAs with 0s, they will be dealt with by titer types
logtiter_layers <- Racmacs:::logtiterTableLayers(merged_map)
logtiter_layers <- lapply(logtiter_layers, \(logtiters) {
  logtiters[is.na(logtiters)] <- 0
  logtiters
})

# Assemble the data list for the model
data_list <- list(
  only_prior = 0,
  N_datasets = numLayers(merged_map),
  N_ags = numAntigens(merged_map),
  N_srs = numSera(merged_map),
  N_sr_groups = numSeraGroups(merged_map),
  sr_groups = as.numeric(srGroups(merged_map)),
  logtiters = simplify2array(logtiter_layers),
  titertypes = simplify2array(Racmacs:::titertypesTableLayers(merged_map))
)

# Decide on starting conditions for the model
sr_group_gmts <- srGroupGMTs(merged_map, ci_method = 'quap')
sr_group_gmts[is.na(sr_group_gmts)] <- 7

data_init <- list(
  sr_group_gmts = sr_group_gmts,
  sr_effects = rep(0, numSera(merged_map)),
  logtiter_error_sigma = rep(0.8, numLayers(merged_map))
)

# Fetch the model
mod <- cmdstan_model(stan_file = 'code/stan_models/titer_magnitude.stan')

# Run the model
mod_optimized <- mod$optimize(
  data = data_list,
  init = list(data_init)
)

# Get the fits
extract_indexnum <- function(variable) as.integer(gsub('.*\\[([0-9]+)\\]$', '\\1', variable))

# Parse the model results
model_result <- list(
  sr_effects = mutate(
    mod_optimized$summary('sr_effects'),
    sr_num = extract_indexnum(variable),
    sr_name = factor(sr_num, labels = srNames(merged_map))
  ),
  dataset_effects = mutate(
    mod_optimized$summary('dataset_effects'),
    dataset_num = extract_indexnum(variable),
    dataset = factor(dataset_num, labels = layerNames(merged_map))
  )
)

# Generate a dataframe with all datasets and adapted titers
# Split out maps
split_maps <- Racmacs:::splitTiterLayers(merged_map)

# Get the map info
long_map_list_info(
  split_maps
) %>%
  mutate(
    ag_name = factor(ag_name, agNames(merged_map)),
    map = factor(map, layerNames(merged_map))
  ) -> map_info

# Format effects
model_result$sr_effects %>%
  rename(sr_effect = estimate) %>%
  select(
    sr_effect, sr_name
  ) -> sr_effects

# Join in effects
map_info %>% left_join(
  sr_effects,
  by = 'sr_name'
) -> map_info

# Add additional columns for the dataset_effects
model_result$dataset_effects %>%
  rename(dataset_effects = estimate) %>%
  select(
    dataset_effects, dataset
  ) %>%
  mutate(
    map = dataset
  ) -> dataset_effects

# Join in effects
map_info %>% left_join(
  dataset_effects,
  by = 'map'
) -> map_info

# For each dataset remove rows from sera and variants that have not been titrated
map_info %>%
  subset(
    titertype != -1
  ) -> map_info_subset

# Save results
saveRDS(map_info, 'data/titer_analyses/titer_magnitude/dataset_magnitude_effect/dataset_magnitude_effect_map_info.rds')
saveRDS(map_info_subset, 'data/titer_analyses/titer_magnitude/dataset_magnitude_effect/dataset_magnitude_effect_map_info_subset.rds')
