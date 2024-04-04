
# Setup workspace
rm(list = ls())
set.seed(100)
library(tidyverse)
library(Racmacs)
library(patchwork)

# Load functions
source('code/metadata/homologous_ags.R')

# Read the map
map <- read.acmap('data/individual_datasets/innsbruck/maps/map-OmicronI+II+III-thresholded-single_exposure-P1m1-no-nds.ace')

sr_groups <- c(
  'B.1.617.2 convalescent', 'B.1.617.2 convalescent', 'B.1.617.2 convalescent', 'B.1.617.2 convalescent', 'B.1.617.2 convalescent', 'B.1.1.7 convalescent', 'B.1.1.7 convalescent', 'B.1.1.7 convalescent', 'B.1.1.7 convalescent', 'B.1.1.7 convalescent', 'B.1.1.7 convalescent', 'B.1.1.7 convalescent', 'B.1.1.7 convalescent', 'B.1.1.7 convalescent', 'B.1.351 convalescent', 'B.1.351 convalescent', 'B.1.351 convalescent', 'B.1.351 convalescent', 'B.1.351 convalescent', 'B.1.351 convalescent', 'mRNA-1273',  'mRNA-1273', 'mRNA-1273', 'mRNA-1273', 'mRNA-1273', 'mRNA-1273', 'mRNA-1273', 'mRNA-1273', 'mRNA-1273', 'mRNA-1273', 'AstraZeneca', 'AstraZeneca', 'AstraZeneca', 'AstraZeneca', 'AstraZeneca', 'AstraZeneca', 'AstraZeneca', 'AstraZeneca', 'AstraZeneca', 'AstraZeneca', 'AstraZeneca-Pfizer', 'AstraZeneca-Pfizer', 'AstraZeneca-Pfizer', 'AstraZeneca-Pfizer', 'AstraZeneca-Pfizer', 'AstraZeneca-Pfizer', 'AstraZeneca-Pfizer', 'AstraZeneca-Pfizer', 'AstraZeneca-Pfizer', 'AstraZeneca-Pfizer', 'Pfizer', 'Pfizer', 'Pfizer', 'Pfizer', 'Pfizer', 'Pfizer', 'Pfizer', 'Pfizer', 'Pfizer', 'Pfizer', 'Pfizer', 'BA.1 convalescent', 'BA.1 convalescent', 'BA.1 convalescent', 'BA.1 convalescent', 'BA.1 convalescent', 'BA.1 convalescent', 'BA.1 convalescent', 'BA.1 convalescent', 'BA.2 convalescent', 'BA.2 convalescent', 'BA.2 convalescent', 'BA.2 convalescent', 'BA.2 convalescent', 'BA.2 convalescent', 'BA.2 convalescent', 'BA.2 convalescent', 'BA.2 convalescent', 'D614G convalescent', 'D614G convalescent', 'D614G convalescent', 'D614G convalescent', 'D614G convalescent', 'D614G convalescent', 'D614G convalescent', 'D614G convalescent', 'D614G convalescent', 'D614G convalescent'
)
srGroups(map) <- factor(sr_groups, levels = unique(sr_groups))


# Exclude certain antigens from higher than homologous comparison
excluded_comp_ags <- c('P.1.1', 'BA.2')

# Set the homologous antigens
srHomologousAgs(map) <- as.list(
  match(homologous_ags_innsbruck[as.character(srGroups(map))], agNames(map))
)

# Determine those individuals with more than 2-fold higher than homologous titers
sr_homologous_logtiters <- unlist(srHomologousLogTiters(map))
sr_max_logtiters <- unname(apply(logtiterTable(map)[!agNames(map) %in% excluded_comp_ags,], 2, max, na.rm = TRUE))
excluded_sr_higher_than_homologous <- srNames(map)[sr_max_logtiters > sr_homologous_logtiters + 1]

# Remove NaNs for sera that have not been titrated against the homologous variants
excluded_sr_higher_than_homologous <- excluded_sr_higher_than_homologous[!is.na(excluded_sr_higher_than_homologous)]

# Set manually judged possible second infections
omicron_second <- c(
  'G776_BA.2 conv._BA.2_NA_NA'
)

other_second <- c()

sera_not_titrated_against_homologous <- c()

manual_excluded_sera <- c(
  omicron_second,
  other_second,
  sera_not_titrated_against_homologous
)

# Save the excluded sera
saveRDS(
  object = list(
    higher_than_homologous = excluded_sr_higher_than_homologous,
    manual = manual_excluded_sera
  ),
  file = 'data/individual_datasets/innsbruck/outliers/outliers.rds'
)

# Create a map with outlier data indicated
possible_second <- c(excluded_sr_higher_than_homologous, manual_excluded_sera)
map_full_no_outliers <- subsetMap(
  map = map,
  sera = !srNames(map) %in% possible_second
)

map_full_no_outliers <- optimizeMap(map_full_no_outliers, 2, 1000)
map_full_no_outliers <- realignMap(map_full_no_outliers, map)

save.acmap(map_full_no_outliers, 'data/individual_datasets/innsbruck/maps/map-OmicronI+II+III-thresholded-single_exposure-P1m1-no-nds-no-outliers.ace')

