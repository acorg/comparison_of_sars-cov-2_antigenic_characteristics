
# Setup workspace
rm(list = ls())
set.seed(100)
library(Racmacs)

mapColors <- read.csv(file = 'data/metadata/map-colors.csv', row.names = 'Antigen', header = TRUE)
mapColorsVec <- unlist(mapColors)
names(mapColorsVec) <- rownames(mapColors)

# Make map from the first set of data
titer_table <- read.titerTable('data/individual_datasets/madison_unpooled/rawdata/individual-hamsters.csv')

ind <- acmap(titer_table = titer_table)
dilutionStepsize(ind) <- 1

ind <- optimizeMap(ind, 2, 500)

pooled <- read.acmap('data/individual_datasets/madison_pooled/maps/map.ace')

ind <- realignMap(ind, pooled)

agSize(ind) <- 18
srSize(ind) <- 10
srOutlineWidth(ind) <- 1

srOutline(ind) <- c(rep('#e7ba52', 6), rep('#393b79', 6))

srGroups(ind) <- factor(c(rep('B.1.351 convalescent', 6), rep('D614G convalescent', 6)))


# Make map from second set of data

# Read titers
titer_table <- readxl::read_excel('data/individual_datasets/madison_unpooled/rawdata/individual hamster serum against 10x virus panel_141221.xlsx', skip = 1)
titer_table <- as.matrix(titer_table)
rownames(titer_table) <- titer_table[,1]
titer_table <- titer_table[,-1]
titer_table <- trimws(titer_table)
titer_table <- t(titer_table)

# Create map
map_141221 <- acmap(
  ag_names = rownames(titer_table),
  sr_names = colnames(titer_table),
  titer_table = titer_table
)

# Standardise the antigen names
agNames(map_141221)[agNames(map_141221) == 'S-614G'] <- 'D614G'
agNames(map_141221)[agNames(map_141221) == 'B.1.526 E484K'] <- 'B.1.526+E484K'
agNames(map_141221)[agNames(map_141221) == 'B.1.526 S477N'] <- 'B.1.526+S477N'
agNames(map_141221)[agNames(map_141221) == 'B.617.2'] <- 'B.1.617.2'

# Standardise the sera names
srNames(map_141221) <- gsub(' #', '-', srNames(map_141221))
srNames(map_141221) <- gsub('Delta', 'B.1.617.2', srNames(map_141221))
srNames(map_141221) <- gsub('Mu', 'B.1.621', srNames(map_141221))
srNames(map_141221) <- gsub('Beta', 'B.1.351', srNames(map_141221))
srNames(map_141221) <- gsub('S-614G', 'D614G', srNames(map_141221))

# Get serum groups
sr_groups <- gsub('(^.*)-.*?$', '\\1 convalescent', srNames(map_141221))
sr_groups <- factor(sr_groups, unique(sr_groups))
srGroups(map_141221) <- sr_groups

# Merge the map with the old data
ind_merged_141221 <- mergeMaps(
  maps = list(
    `individual-hamsters.csv` = ind,
    `individual hamster serum against 10x virus panel_141221.xlsx` = map_141221
  )
)

# Optimize the map
ind_merged_141221 <- optimizeMap(
  map = ind_merged_141221,
  number_of_dimensions = 2,
  number_of_optimizations = 500,
  minimum_column_basis = 'none'
)

# Set map styles
agSize(ind_merged_141221) <- 18
srSize(ind_merged_141221) <- 10
srOutlineWidth(ind_merged_141221) <- 1
agFill(ind_merged_141221) <- mapColorsVec[agNames(ind_merged_141221)]
srOutline(ind_merged_141221) <- mapColorsVec[as.character(srGroups(ind_merged_141221))]
ptDrawingOrder(ind_merged_141221) <- rev(ptDrawingOrder(ind_merged_141221))

ind_merged_141221 <- realignMap(ind_merged_141221, pooled)
ind_merged_141221 <- reflectMap(ind_merged_141221, axis = 'y')
ind_merged_141221 <- rotateMap(ind_merged_141221, 50)

# Save the map
save.acmap(ind_merged_141221, 'data/individual_datasets/madison_unpooled/maps/map.ace')

