
# Setup workspace
rm(list = ls())
set.seed(100)
library(tidyverse)
library(Racmacs)

library(titertools)

# Read in the map
map_full_no_outliers <- read.acmap('data/individual_datasets/fda/maps/map_no_outliers.ace')

# Calculate individual effects
# Get serum group data
individual_effects <- tibble(
  sr_name = character(0),
  sr_individual_effect = numeric(0)
)

ag_means_individual_effects <- tibble(
  sr_group = factor(NULL, levels(srGroups(map_full_no_outliers))),
  ag_name = character(0),
  ag_gmt_accounting_for_individual_effect = numeric(0)
)

# Calculate sr group individual effects
for (sr_group in unique(srGroups(map_full_no_outliers))) {

  message(sr_group)

  sr_group_sr_names <- srNames(map_full_no_outliers)[srGroups(map_full_no_outliers) == sr_group]
  titers <- titerTable(map_full_no_outliers)[,srGroups(map_full_no_outliers) == sr_group, drop = F]

  # Calculate the fit
  individual_effects_result <- titertools::estimate_sr_effects(
    titers = titers,
    dilution_stepsize = 0,
    sigma = 0.62,
    ci_method = 'quap'
  )

  # Construct tables
  sr_group_individual_effects <- tibble(
    sr_name = sr_group_sr_names,
    sr_individual_effect = individual_effects_result[sprintf('sr_effects[%s]', seq_len(ncol(titers))), 'estimate']
  )

  sr_group_ag_means_individual_effects <- tibble(
    ag_name = agNames(map_full_no_outliers),
    sr_group = factor(sr_group, levels(srGroups(map_full_no_outliers))),
    ag_gmt_accounting_for_individual_effect = individual_effects_result[sprintf('ag_means[%s]', seq_len(nrow(titers))), 'estimate']
  )

  individual_effects <- bind_rows(
    individual_effects,
    sr_group_individual_effects
  )

  ag_means_individual_effects <- bind_rows(
    ag_means_individual_effects,
    sr_group_ag_means_individual_effects
  )

}

# Save the individual effects
saveRDS(individual_effects, file = 'data/individual_datasets/fda/individual_effects/individual_effects.rds')
saveRDS(ag_means_individual_effects, file = 'data/individual_datasets/fda/individual_effects/ag_gmt_individual_effects.rds')
