
# Setup workspace
rm(list = ls())
set.seed(100)

library(Racmacs)

map <- read.acmap('data/individual_datasets/innsbruck/maps/map-OmicronI+II+III-thresholded-single_exposure-P1m1.ace')

# Remove sera with less than 3 detectable titers
ndetectable <- colSums(titerTable(map) != '*' & !grepl('<', titerTable(map)))
sera_subset <- srNames(map)[ndetectable >= 3]

map_no_nd <- subsetMap(
  map, 
  sera = sera_subset
)

map_no_nd <- optimizeMap(map_no_nd, 2, 1000)

map_no_nd <- realignMap(map_no_nd, map)

save.acmap(map_no_nd, 'data/individual_datasets/innsbruck/maps/map-OmicronI+II+III-thresholded-single_exposure-P1m1-no-nds.ace')