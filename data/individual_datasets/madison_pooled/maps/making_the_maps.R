
# Setup workspace
rm(list = ls())
set.seed(100)
library(Racmacs)

mapColors <- read.csv(file = 'data/metadata/map-colors.csv', row.names = 'Antigen', header = TRUE)
mapColorsVec <- unlist(mapColors)
names(mapColorsVec) <- rownames(mapColors)
nameEquivalents <- read.csv(file = 'data/metadata/name-equivalence.csv', row.names = 'origName', header = TRUE)

# Make map from the first set of data received
m <- read.acmap('data/individual_datasets/duke/maps/map_no_outliers.ace')
titer_table <- read.titerTable('data/individual_datasets/madison_pooled/rawdata/madison-hamster.csv')

firstMap <- acmap(titer_table = titer_table)
dilutionStepsize(firstMap) <- 1
firstMap <- optimizeMap(firstMap, 2, 500)

firstMap <- realignMap(firstMap, m)

agSize(firstMap) <- 18
srSize(firstMap) <- 10
srOutlineWidth(firstMap) <- 1

agFill(firstMap)[agNames(firstMap) == 'B.1.1.7+E484K'] <- mapColors['B.1.1.7+E484K',]
agFill(firstMap)[agNames(firstMap) == 'B.1.526+E484K'] <- mapColors['B.1.526+E484K',]
agFill(firstMap)[agNames(firstMap) == 'B.1.526+S477N'] <- mapColors['B.1.526+S477N',]
agFill(firstMap)[agNames(firstMap) == 'C.36.3'] <- mapColors['C.36.3',]
agFill(firstMap)[agNames(firstMap) == 'B.1.617.1'] <- mapColors['B.1.617.1',]
agFill(firstMap)[agNames(firstMap) == 'B.1.617.2'] <- mapColors['B.1.617.2',]
agFill(firstMap)[agNames(firstMap) == 'D614G'] <- mapColors['D614G',]
agFill(firstMap)[agNames(firstMap) == 'P.1'] <- mapColors['P.1',]
agFill(firstMap)[agNames(firstMap) == 'B.1.351'] <- mapColors['B.1.351',]
agFill(firstMap)[agNames(firstMap) == 'B.1.1.7'] <- mapColors['B.1.1.7',]

srOutline(firstMap) <- unlist(lapply(srNames(firstMap), function(i) mapColors[i,]), use.names=FALSE)

srGroups(firstMap) <- factor(c('B.1.617.1 convalescent', 'B.1.617.2 convalescent', 'D614G convalescent', 'P.1 convalescent', 'B.1.351 convalescent', 'B.1.1.7 convalescent'))

# Make map from the second set of data received on 2021-08-31
titer_table <- read.titerTable('data/individual_datasets/madison_pooled/rawdata/madison-hamster-update-210831.csv')

kaw2 <- acmap(titer_table = titer_table)
dilutionStepsize(kaw2) <- 1
kaw2 <- optimizeMap(kaw2, 2, 500)

kaw2 <- applyPlotspec(kaw2, firstMap)
kaw2 <- realignMap(kaw2, firstMap)

agSize(kaw2) <- 18
srSize(kaw2) <- 10
srOutlineWidth(kaw2) <- 1

agFill(kaw2)[agNames(kaw2) == 'B.1.621'] <- mapColors['B.1.621',]
agFill(kaw2)[agNames(kaw2) == '614D'] <- mapColors['614D',]

srOutline(kaw2) <- unlist(lapply(srNames(kaw2), function(i) mapColors[i,]), use.names=FALSE)

srGroups(kaw2) <- factor(c('B.1.617.1 convalescent', 'C.36.3 convalescent', 'B.1.351 convalescent', 'B.1.621 convalescent', 'B.1.617.2 convalescent', '614D convalescent'))

ag_groups <- unlist(lapply(agNames(kaw2), function(i) nameEquivalents[i,]), use.names=FALSE)
agGroups(kaw2) <- factor(ag_groups)


# Make a map from the merge of the two datasets
kawMerge <- mergeMaps(list(firstMap, kaw2), method = 'reoptimized-merge', number_of_dimensions = 2, number_of_optimizations = 500,
                      minimum_column_basis = 'none')

kawMerge <- realignMap(kawMerge, firstMap)

srGroups(kawMerge) <- factor(c('B.1.617.1 convalescent', 'B.1.617.2 convalescent', '614D convalescent', 'P.1 convalescent', 'B.1.351 convalescent', 'B.1.1.7 convalescent', 'C.36.3 convalescent', 'B.1.621 convalescent', '614D convalescent'))

ag_groups <- unlist(lapply(agNames(kawMerge), function(i) nameEquivalents[i,]), use.names=FALSE)
agGroups(kawMerge) <- factor(ag_groups)

# Make the reactivity adjusted map
# Adjust antigen reactivity
agReactivityAdjustments(kawMerge) <- c(0, 0, 0, -3, 0, -3, 0, 0, 0, 0, 0, 0)

# Optimize the map
pooled <- optimizeMap(
  kawMerge,
  number_of_dimensions = 2,
  number_of_optimizations = 1000,
  minimum_column_basis = 'none'
)

# Set antigen colors
agFill(pooled)[agNames(pooled) == 'B.1.1.7+E484K'] <- '#637939'

# Set styles
srOutlineWidth(pooled) <- 1
srSize(pooled) <- 10
agSize(pooled) <- 18

smaller_ags <- c(
  'B.1.1.7+E484K'
)
agSize(pooled)[agNames(pooled) %in% smaller_ags] <- 12

# Set map orientation
reference_map <- read.acmap('data/individual_datasets/duke/maps/map_no_outliers.ace')
pooled <- realignMap(pooled, reference_map)

# Name the map
mapName(pooled) <- 'madison_pooled'

# Save the map
save.acmap(pooled, 'data/individual_datasets/madison_pooled/maps/map.ace')
