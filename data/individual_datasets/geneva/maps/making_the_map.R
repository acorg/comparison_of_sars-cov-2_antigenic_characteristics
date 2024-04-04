
# Setup workspace
rm(list = ls())
set.seed(100)
library(Racmacs)

# Read in the titer data
titer_table <- read.titerTable('data/individual_datasets/geneva/rawdata/figure3-data-reformatted.csv')

map <- acmap(titer_table = titer_table)
dilutionStepsize(map) <- 0

map <- optimizeMap(map, 2, 1000)

agSize(map) <- 18
srSize(map) <- 10
srOutlineWidth(map) <- 1

ptDrawingOrder(map) <- rev(ptDrawingOrder(map))

agNames(map) <- c('D614G', 'B.1.1.7', 'B.1.351', 'P.1', 'B.1.617.2', 'P.2', 'BA.1')
agFill(map) <- c('#393b79', '#637939', '#e7ba52', '#7b4173', '#d18652', '#be7cf8', '#EF3737')
sr_groups <- c('D614G convalescent', 'D614G convalescent', 'D614G convalescent', 'D614G convalescent', 
               'D614G convalescent', 'D614G convalescent', 'D614G convalescent', 'D614G convalescent', 
               'D614G convalescent', 'D614G convalescent', 'D614G convalescent', 'D614G convalescent', 
               'D614G convalescent', 'D614G convalescent', 'D614G convalescent', 'D614G convalescent', 
               'D614G convalescent', 'D614G convalescent', 'D614G convalescent', 'D614G convalescent', 
               'D614G convalescent', 'D614G convalescent', 'D614G convalescent', 'D614G convalescent', 
               'D614G convalescent', 'D614G convalescent', 'D614G convalescent', 'D614G convalescent', 
               'D614G convalescent', 'D614G convalescent', 'D614G convalescent', 'D614G convalescent', 
               'D614G convalescent', 'D614G convalescent', 'B.1.1.7 convalescent', 'B.1.1.7 convalescent', 
               'B.1.1.7 convalescent', 'B.1.1.7 convalescent', 'B.1.1.7 convalescent', 'B.1.1.7 convalescent', 
               'B.1.1.7 convalescent', 'B.1.1.7 convalescent', 'B.1.1.7 convalescent', 'B.1.1.7 convalescent', 
               'B.1.1.7 convalescent', 'B.1.1.7 convalescent', 'B.1.351 convalescent', 'B.1.351 convalescent', 
               'B.1.351 convalescent', 'B.1.351 convalescent', 'B.1.351 convalescent', 'B.1.351 convalescent', 
               'B.1.351 convalescent', 'B.1.351 convalescent', 'P.1 convalescent', 'P.1 convalescent', 
               'P.1 convalescent', 'P.1 convalescent', 'P.1 convalescent', 'P.1 convalescent', 
               'P.1 convalescent', 'P.1 convalescent', 'P.1 convalescent', 'P.1 convalescent', 
               'B.1.617.2 convalescent', 'B.1.617.2 convalescent', 'B.1.617.2 convalescent', 'B.1.617.2 convalescent', 
               'B.1.617.2 convalescent', 'B.1.617.2 convalescent', 'B.1.617.2 convalescent', 'B.1.617.2 convalescent', 
               'B.1.617.2 convalescent', 'B.1.617.2 convalescent', 'mRNA-1273', 'mRNA-1273', 'mRNA-1273', 'mRNA-1273', 
               'mRNA-1273', 'mRNA-1273', 'mRNA-1273', 'mRNA-1273', 'mRNA-1273', 'mRNA-1273', 'mRNA-1273', 'mRNA-1273', 
               'mRNA-1273', 'mRNA-1273', 'mRNA-1273', 'mRNA-1273')
srGroups(map) <- factor(sr_groups, levels = c('mRNA-1273', 'D614G convalescent', 'B.1.1.7 convalescent',
                                              'B.1.351 convalescent', 'P.1 convalescent',
                                              'B.1.617.2 convalescent'))
srOutline(map)[srGroups(map) == 'D614G convalescent'] <- '#333333'
srOutline(map)[srGroups(map) == 'B.1.1.7 convalescent'] <- '#637939'
srOutline(map)[srGroups(map) == 'B.1.351 convalescent'] <- '#e7ba52'
srOutline(map)[srGroups(map) == 'B.1.617.2 convalescent'] <- '#d18652'
srOutline(map)[srGroups(map) == 'P.1 convalescent'] <- '#7b4173'
srOutline(map)[srGroups(map) == 'mRNA-1273'] <- 'grey'

map <- reflectMap(map, axis = 'x')
map <- rotateMap(map, 0.5)

map <- translateMap(map, c(1, 1.2))

save.acmap(map, 'data/individual_datasets/geneva/maps/map.ace')

