
# Setup workspace
rm(list = ls())
set.seed(100)
library(Racmacs)

titer_table <- read.titerTable('data/individual_datasets/amc/rawdata/titer-table.csv')

# Sera with less than 3 detectable titers were removed. 
# Sera with unknown infection history (BA.1 / BA.2 were removed).
# Sera where homologous antigen wasn't titrated were removed (Alpha+E484K).
# Beta without the deletion was removed.

map <- acmap(titer_table = titer_table)
dilutionStepsize(map) <- 0

agSize(map) <- 18
srSize(map) <- 10
srOutlineWidth(map) <- 1

ptDrawingOrder(map) <- rev(ptDrawingOrder(map))

agNames(map) <- c('D614G', 'B.1.1.7', 'B.1.351-subst', 'B.1.351', 'P.1', 'B.1.617.2', 'BA.1', 'BA.2')
agFill(map) <- c('#393b79', '#637939', '#e7ba52', '#e7ba52', '#7b4173', '#d18652', '#EF3737', '#d10fa2')

# Remove titrations against the Beta variant with the L242H, R246I substitution
# Remove titrations against the sera with unknown BA.1 or BA.2 first infection
# Remove titrations against sera where the homologous variant wasn't titrated (B.1.1.7+E484K)
map_subset <- subsetMap(map, antigens = agNames(map)[agNames(map) != 'B.1.351-subst'], sera = srNames(map)[!srNames(map) %in% c('COSCA-352', 'COSCA-353', 'COSCA-316', 'COSCA-320')])

sr_groups <- c('D614G convalescent', 'D614G convalescent', 'D614G convalescent', 'D614G convalescent', 'D614G convalescent', 'D614G convalescent', 'D614G convalescent', 'D614G convalescent', 'D614G convalescent', 'D614G convalescent', 'D614G convalescent', 'D614G convalescent', 'D614G convalescent', 'D614G convalescent', 'D614G convalescent', 'D614G convalescent', 'D614G convalescent', 'D614G convalescent', 'D614G convalescent', 'D614G convalescent', 'B.1.1.7 convalescent', 'B.1.1.7 convalescent', 'B.1.1.7 convalescent', 'B.1.1.7 convalescent', 'B.1.1.7 convalescent', 'B.1.1.7 convalescent', 'B.1.1.7 convalescent', 'B.1.1.7 convalescent', 'B.1.1.7 convalescent', 'B.1.351 convalescent', 'B.1.351 convalescent', 'B.1.351 convalescent', 'P.1 convalescent', 'P.1 convalescent', 'P.1 convalescent', 'P.1 convalescent', 'B.1.617.2 convalescent', 'B.1.617.2 convalescent', 'B.1.617.2 convalescent', 'B.1.617.2 convalescent', 'B.1.617.2 convalescent', 'B.1.617.2 convalescent', 'B.1.617.2 convalescent', 'B.1.617.2 convalescent', 'B.1.617.2 convalescent', 'B.1.617.2 convalescent', 'B.1.617.2 convalescent', 'B.1.351 convalescent', 'B.1.351 convalescent', 'B.1.351 convalescent', 'B.1.351 convalescent', 'B.1.351 convalescent', 'BA.1 convalescent', 'BA.1 convalescent', 'BA.1 convalescent', 'BA.1 convalescent', 'BA.1 convalescent', 'BA.2 convalescent', 'BA.2 convalescent', 'BA.2 convalescent', 'BA.1 convalescent', 'BA.2 convalescent')

srGroups(map_subset) <- factor(sr_groups, levels = c('D614G convalescent', 'B.1.1.7 convalescent', 'B.1.351 convalescent',
                                                     'P.1 convalescent', 'B.1.617.2 convalescent', 'BA.1 convalescent',
                                                     'BA.2 convalescent'))

srOutline(map_subset)[srGroups(map_subset) == 'D614G convalescent'] <- '#333333'
srOutline(map_subset)[srGroups(map_subset) == 'B.1.1.7 convalescent'] <- '#637939'
srOutline(map_subset)[srGroups(map_subset) == 'B.1.351 convalescent'] <- '#e7ba52'
srOutline(map_subset)[srGroups(map_subset) == 'P.1 convalescent'] <- '#7b4173'
srOutline(map_subset)[srGroups(map_subset) == 'B.1.617.2 convalescent'] <- '#d18652'
srOutline(map_subset)[srGroups(map_subset) == 'BA.1 convalescent'] <- '#EF3737'
srOutline(map_subset)[srGroups(map_subset) == 'BA.2 convalescent'] <- '#d10fa2'

# Remove sera with less than 3 detectable titers
ndetectable <- colSums(titerTable(map_subset) != '*' & !grepl('<', titerTable(map_subset)))
sera_subset <- srNames(map_subset)[ndetectable >= 3]

map_subset_no_nd <- subsetMap(
  map_subset, 
  sera = sera_subset
)

map_subset_no_nd <- optimizeMap(map_subset_no_nd, 2, 1000)

map <- rotateMap(map_subset_no_nd, -26)


save.acmap(map, 'data/individual_datasets/amc/maps/map.ace')
