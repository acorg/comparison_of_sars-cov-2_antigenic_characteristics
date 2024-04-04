
# Setup workspace
rm(list = ls())
library(tidyverse)
library(Racmacs)
set.seed(100)

# Load the map for comparisons
duke <- read.acmap('data/individual_datasets/duke/maps/map_no_outliers.ace')

# Load the new omicron data
titers <- read.titerTable('data/individual_datasets/madison_frnt/rawdata/220410-Hamster-serum-FRNT.csv')
sr_groups_info <- read.csv('data/individual_datasets/madison_frnt/rawdata/sr_groups.csv', stringsAsFactors = FALSE)

# Create a map from all the data
frnt <- acmap(
  ag_names = rownames(titers),
  sr_names = colnames(titers),
  titer_table = titers
)

dilutionStepsize(frnt) <- 0

sr_groups <- sr_groups_info$Serum.group[match(srNames(frnt), sr_groups_info$Serum)]
srGroups(frnt) <- factor(
  sr_groups,
  levels = c(
    '614D convalescent',
    'B.1.617.2 convalescent',
    'BA.1 convalescent',
    'BA.1.1 convalescent',
    'BA.2 convalescent'
  )
)

# Set antigen colors
ag_colors_info <- read.csv('data/metadata/ag_colors.csv', stringsAsFactors = FALSE)
ag_colors <- ag_colors_info$Color[match(agNames(frnt), ag_colors_info$Antigen)]
agFill(frnt) <- ag_colors

# Set serum colors
sr_colors_info <- read.csv('data/metadata/sr_colors.csv', stringsAsFactors = FALSE)
sr_colors <- sr_colors_info$Color[match(srGroups(frnt), sr_colors_info$Serum)]
srOutline(frnt) <- sr_colors

# Set styles
srOutlineWidth(frnt) <- 1
srSize(frnt) <- 10
agSize(frnt) <- 18

agSize(frnt)[agNames(frnt) == 'BA.1.1'] <- 12

# Remove sera and antigens with less than 3 detectable titers
ndetectable_sr <- colSums(titerTable(frnt) != '*' & !grepl('<', titerTable(frnt)))
ndetectable_ag <- rowSums(titerTable(frnt) != '*' & !grepl('<', titerTable(frnt)))

frnt_nd_subset <- subsetMap(frnt, sera = ndetectable_sr > 2)
frnt_nd_subset <- subsetMap(frnt_nd_subset, antigens = ndetectable_ag > 2)

# Optimize the map
frnt_nd_subset <- optimizeMap(
  frnt_nd_subset,
  number_of_dimensions = 2,
  number_of_optimizations = 1000,
  minimum_column_basis = 'none'
)

# Realign maps
frnt_nd_subset <- realignMap(frnt_nd_subset, duke)
frnt_nd_subset <- translateMap(frnt_nd_subset, c(-0.5, 1))

# Name the map
mapName(frnt_nd_subset) <- 'madison_frnt'

# Set point drawing order
ptDrawingOrder(frnt_nd_subset) <- rev(ptDrawingOrder(frnt_nd_subset))

# Save the maps
save.acmap(frnt_nd_subset, 'data/individual_datasets/madison_frnt/maps/map.ace')
