
# Setup workspace
rm(list = ls())
library(tidyverse)
library(Racmacs)
set.seed(100)

# Read in the raw titers and serum group info
titers <- read.titerTable('data/individual_datasets/mt_sinai/rawdata/table-human.csv')
sr_groups_info <- read.csv('data/individual_datasets/mt_sinai/rawdata/sr_groups.csv', stringsAsFactors = FALSE)

# Create a map from all the data
map_full <- acmap(
  ag_names = rownames(titers),
  sr_names = colnames(titers),
  titer_table = titers
)
dilutionStepsize(map_full) <- 0

# Set serum groups
sr_groups <- sr_groups_info$Serum.group[match(srNames(map_full), sr_groups_info$Serum)]
srGroups(map_full) <- factor(
  sr_groups,
  levels = c(
    'mRNA-1273',
    'B.1.351 convalescent'
  )
)

# Edit ag names
agNames(map_full)[agNames(map_full) == 'B.1.526 (E484K)'] <- 'B.1.526+E484K'
agNames(map_full)[agNames(map_full) == 'B.1.526 (S477)'] <- 'B.1.526+S477N'
agNames(map_full)[agNames(map_full) == 'WA1'] <- '614D'

# Set antigen colors
ag_colors_info <- read.csv('data/metadata/ag_colors.csv', stringsAsFactors = FALSE)
ag_colors <- ag_colors_info$Color[match(agNames(map_full), ag_colors_info$Antigen)]
agFill(map_full) <- ag_colors

# Set serum colors
sr_colors_info <- read.csv('data/metadata/sr_colors.csv', stringsAsFactors = FALSE)
sr_colors <- sr_colors_info$Color[match(srGroups(map_full), sr_colors_info$Serum)]
srOutline(map_full) <- sr_colors

# Set styles
srOutlineWidth(map_full) <- 1
srSize(map_full) <- 10
agSize(map_full) <- 18

smaller_ags <- c(
  'B.1.1.7+E484K'
)
agSize(map_full)[agNames(map_full) %in% smaller_ags] <- 12
ptDrawingOrder(map_full) <- rev(ptDrawingOrder(map_full))

# Subset the map to sera with >= 3 detectable titers
ndetectable <- colSums(titerTable(map_full) != '*' & !grepl('<', titerTable(map_full)))
map_subset <- subsetMap(map_full, sera = ndetectable > 2)

# Optimize the map
map_subset <- optimizeMap(
  map_subset,
  number_of_dimensions = 2,
  number_of_optimizations = 1000,
  minimum_column_basis = 'none'
)

# Set map orientation
reference_map <- read.acmap('data/individual_datasets/duke/maps/map_no_outliers.ace')
map_subset <- realignMap(map_subset, reference_map)
map_subset <- reflectMap(map_subset, axis = 'y')
map_subset <- rotateMap(map_subset, 45)
map_subset <- translateMap(map_subset, c(3.25, -1.4))

# Name the map
mapName(map_subset) <- 'mt_sinai'

# Save the map
save.acmap(map_subset, 'data/individual_datasets/mt_sinai/maps/map.ace')

