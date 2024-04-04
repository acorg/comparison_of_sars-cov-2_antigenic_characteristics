
# Setup workspace
rm(list = ls())
set.seed(100)
library(Racmacs)

# PRNT map
# Read in the raw titers and serum group info
titers <- read.titerTable('data/individual_datasets/emc/rawdata/table-220328.csv')
sr_groups_info <- read.csv('data/individual_datasets/emc/rawdata/sr_groups-220328.csv', stringsAsFactors = FALSE)

# Add serum group information
titer_serum_groups <- sr_groups_info$Serum.group[match(colnames(titers), sr_groups_info$Serum)]

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
    'D614G convalescent',
    'B.1.1.7 convalescent',
    'B.1.351 convalescent',
    'P.1 convalescent',
    'P.2 convalescent',
    'B.1.617.2 convalescent',
    'B.1.621 convalescent',
    'BA.1 convalescent'
  )
)

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
  'AY.4.2'
)
agSize(map_full)[agNames(map_full) %in% smaller_ags] <- 12
ptDrawingOrder(map_full) <- rev(ptDrawingOrder(map_full))

# Subset the map to sera with >= 3 detectable titers
ndetectable <- colSums(titerTable(map_full) != '*' & !grepl('<', titerTable(map_full)))
map_subset <- subsetMap(map_full, sera = ndetectable > 2)

# Remove P.1 sera
map_subset <- removeSera(map_subset, c('Gamma-1', 'Gamma-2', 'Gamma-3', 'Gamma-4'))

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
map_subset <- rotateMap(map_subset, -54)

# Name the map
mapName(map_subset) <- 'emc'

# Save the map
save.acmap(map_subset, 'data/individual_datasets/emc/maps/emc-prnt.ace')


# Pseudo Vero map
titer_table <- read.titerTable('data/individual_datasets/emc/rawdata/pseudo-vero.csv')

map <- acmap(titer_table = titer_table)
dilutionStepsize(map) <- 0

map <- optimizeMap(map, 2, 1000)

agSize(map) <- 18
srSize(map) <- 10
srOutlineWidth(map) <- 1

ptDrawingOrder(map) <- rev(ptDrawingOrder(map))

agNames(map) <- c('614D', 'D614G', 'B.1.1.7', 'B.1.351', 'B.1.617.2', 'B.1.617.1', 'BA.1')
agFill(map) <- c('#393b79', '#393b79', '#637939', '#e7ba52', '#d18652', '#ad494a', '#EF3737')
sr_groups <- c('D614G convalescent', 'D614G convalescent', 'D614G convalescent',
               'D614G convalescent', 'B.1.1.7 convalescent', 'B.1.1.7 convalescent',
               'B.1.1.7 convalescent', 'B.1.1.7 convalescent', 'B.1.351 convalescent',
               'B.1.351 convalescent',  'B.1.351 convalescent',  'B.1.351 convalescent',
               'B.1.617.2 convalescent', 'B.1.617.2 convalescent', 'B.1.617.2 convalescent',
               'B.1.617.2 convalescent', 'BA.1 convalescent',  'BA.1 convalescent',
               'BA.1 convalescent',  'BA.1 convalescent')
srGroups(map) <- factor(sr_groups, levels = c('D614G convalescent', 'B.1.1.7 convalescent', 'B.1.351 convalescent',
                                              'B.1.617.2 convalescent', 'BA.1 convalescent'))
srOutline(map)[srGroups(map) == 'D614G convalescent'] <- '#333333'
srOutline(map)[srGroups(map) == 'B.1.1.7 convalescent'] <- '#637939'
srOutline(map)[srGroups(map) == 'B.1.351 convalescent'] <- '#e7ba52'
srOutline(map)[srGroups(map) == 'B.1.617.2 convalescent'] <- '#d18652'
srOutline(map)[srGroups(map) == 'BA.1 convalescent'] <- '#EF3737'

map <- rotateMap(map, -1)

save.acmap(map, 'data/individual_datasets/emc/maps/pseudo-vero.ace')



# Pseudo Calu-3 map
titer_table <- read.titerTable('data/individual_datasets/emc/rawdata/pseudo-calu.csv')

map <- acmap(titer_table = titer_table)
dilutionStepsize(map) <- 0

map <- optimizeMap(map, 2, 1000)

agSize(map) <- 18
srSize(map) <- 10
srOutlineWidth(map) <- 1

ptDrawingOrder(map) <- rev(ptDrawingOrder(map))

agNames(map) <- c('614D', 'D614G', 'B.1.1.7', 'B.1.351', 'B.1.617.2', 'B.1.617.1', 'BA.1')
agFill(map) <- c('#393b79', '#393b79', '#637939', '#e7ba52', '#d18652', '#ad494a', '#EF3737')
sr_groups <- c('D614G convalescent', 'D614G convalescent', 'D614G convalescent',
               'D614G convalescent', 'B.1.1.7 convalescent', 'B.1.1.7 convalescent',
               'B.1.1.7 convalescent', 'B.1.1.7 convalescent', 'B.1.351 convalescent',
               'B.1.351 convalescent',  'B.1.351 convalescent',  'B.1.351 convalescent',
               'B.1.617.2 convalescent', 'B.1.617.2 convalescent', 'B.1.617.2 convalescent',
               'B.1.617.2 convalescent', 'BA.1 convalescent',  'BA.1 convalescent',
               'BA.1 convalescent',  'BA.1 convalescent')
srGroups(map) <- factor(sr_groups, levels = c('D614G convalescent', 'B.1.1.7 convalescent', 'B.1.351 convalescent',
                                              'B.1.617.2 convalescent', 'BA.1 convalescent'))
srOutline(map)[srGroups(map) == 'D614G convalescent'] <- '#333333'
srOutline(map)[srGroups(map) == 'B.1.1.7 convalescent'] <- '#637939'
srOutline(map)[srGroups(map) == 'B.1.351 convalescent'] <- '#e7ba52'
srOutline(map)[srGroups(map) == 'B.1.617.2 convalescent'] <- '#d18652'
srOutline(map)[srGroups(map) == 'BA.1 convalescent'] <- '#EF3737'

mV <- read.acmap('data/individual_datasets/emc/maps/pseudo-vero.ace')
map <- realignMap(map, mV)

map <- rotateMap(map, -1)
view(map)

save.acmap(map, 'data/individual_datasets/emc/maps/pseudo-calu.ace')
