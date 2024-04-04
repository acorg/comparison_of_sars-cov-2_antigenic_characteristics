
# Setup workspace
rm(list = ls())
set.seed(100)
library(Racmacs)
library(tidyverse)

mapColors <- read.csv(file = 'data/metadata/map-colors.csv', row.names = 'Antigen', header = TRUE)
nameEquivalents <- read.csv(file = 'data/metadata/name-equivalence.csv', row.names = 'origName', header = TRUE)

# Infer maps for individual tables
titer_table1 <- read.titerTable('data/individual_datasets/emory/rawdata/table-run-1.csv')
map1 <- acmap(titer_table = titer_table1)
dilutionStepsize(map1) <- 0

titer_table2 <- read.titerTable('data/individual_datasets/emory/rawdata/table-run-2.csv')
map2 <- acmap(titer_table = titer_table2)
dilutionStepsize(map2) <- 0

titer_table3 <- read.titerTable('data/individual_datasets/emory/rawdata/table-run-3.csv')
map3 <- acmap(titer_table = titer_table3)
dilutionStepsize(map3) <- 0

titer_table4 <- read.titerTable('data/individual_datasets/emory/rawdata/table-run-4.csv')
map4 <- acmap(titer_table = titer_table4)
dilutionStepsize(map4) <- 0

titer_table5 <- read.titerTable('data/individual_datasets/emory/rawdata/table-run-5.csv')
map5 <- acmap(titer_table = titer_table5)
dilutionStepsize(map5) <- 0

titer_table6 <- read.titerTable('data/individual_datasets/emory/rawdata/table-run-6.csv')
map6 <- acmap(titer_table = titer_table6)
dilutionStepsize(map6) <- 0


# Merge the maps
mergedMap <- mergeMaps(
  list(map1, map2, map3, map4, map5, map6),
  method = 'table',
  number_of_dimensions=2,
  number_of_optimizations='1000',
  minimum_column_basis = 'none'
)

agNames(mergedMap) <- rownames(titer_table1)
srNames(mergedMap) <- colnames(titer_table1)

dilutionStepsize(mergedMap) <- 0

# Remove outliers
titerTable(mergedMap)['B.1.617.2', '867485368'] <- '*'
titerTable(mergedMap)['B.1.1.7_E484K', '762511414'] <- '*'

# Remove sera with <3 detectable titers
mergedMap <- removeSera(mergedMap, c('751630239'))

mergedMap <- optimizeMap(
    map                     = mergedMap,
    number_of_dimensions    = 2,
    number_of_optimizations = 1000,
    minimum_column_basis    = 'none',
    )

ptDrawingOrder(mergedMap) <- rev(seq_len(numPoints(mergedMap)))

ag_groups <- unlist(lapply(agNames(mergedMap), function(i) nameEquivalents[i,]), use.names=FALSE)
agGroups(mergedMap) <- factor(ag_groups)


# Apply styling
srSize(mergedMap) <- 10
srOutlineWidth(mergedMap) <- 1
agSize(mergedMap) <- 18
agSize(mergedMap)[agNames(mergedMap) == 'B.1.1.7_E484K'] <- 12
agSize(mergedMap)[agNames(mergedMap) == 'B.1.617.2_K417N'] <- 12


agFill(mergedMap) <- unlist(lapply(agNames(mergedMap), function(i) mapColors[i,]), use.names=FALSE)

agFill(mergedMap)[agNames(mergedMap) == 'B.1.1.7_E484K'] <- mapColors['B.1.1.7',]
srOutline(mergedMap) <- c(rep(mapColors['WT vac mRNA', ], 17), rep(mapColors['WT convalescent', ], 17), rep(mapColors['B.1.351 convalescent', ], 16), rep(mapColors['B.1.617.2 convalescent', ], 14))

sr_groups <- rep(', numSera(mergedMap))
sr_groups[grepl('#333333', srOutline(mergedMap), fixed = T)] <- 'D614G convalescent'
sr_groups[grepl('#d18652', srOutline(mergedMap), fixed = T)] <- 'B.1.617.2 convalescent'
sr_groups[grepl('#e7ba52', srOutline(mergedMap), fixed = T)] <- 'B.1.351 convalescent'
sr_groups[grepl('grey', srOutline(mergedMap), fixed = T)]    <- 'mRNA-1273'
srGroups(mergedMap) <- factor(sr_groups)

# Adjust antigen reactivity
agReactivityAdjustments(mergedMap) <- c(0, 0, 0, -1, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

mergedMapAdjusted <- optimizeMap(
    map                     = mergedMap,
    number_of_dimensions    = 2,
    number_of_optimizations = 1000,
    minimum_column_basis    = 'none',

    )

save.acmap(mergedMapAdjusted, 'data/individual_datasets/emory/maps/map.ace')
