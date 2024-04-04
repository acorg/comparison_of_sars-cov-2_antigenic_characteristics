
# Setup workspace
rm(list = ls())
set.seed(100)
library(Racmacs)

# Read in the map with additional antigens
map <- read.acmap('data/individual_datasets/charite/maps/map-continuous-fixbottom-90-corrected.ace')

# Subset the map to remove Omicron variants newer than BA.5, and BA.5 convalescent sera.
map_subset <- subsetMap(
  map, 
  sera = srNames(map)[!srNames(map) %in% c('9.1', '9.2', '9.3')],
  antigens = agNames(map)[!agNames(map) %in% c('XBB.2', 'BN.1.3.1', 'BF.7', 'BQ.1.18')]
)

map_subset <- optimizeMap(map_subset, 2, 1000)
map_subset <- realignMap(map_subset, map)
map_subset <- rotateMap(map_subset, -34)

save.acmap(map_subset, 'data/individual_datasets/charite/maps/map.ace')

