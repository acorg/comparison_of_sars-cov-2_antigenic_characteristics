
library(Racmacs)

## Make the map without removing the putative beta sera
map <- read.acmap('data/individual_datasets/fda/rawdata/fifth_packet_primary_sera_ant_map_deidentified.ace')

map <- removeAntigens(map, c('E484K', 'E484K+N501Y', 'K417N', 'L452R', 'N501Y', 'R346K', 'T478K'))
map <- removeSera(map, c('Boca-03-10', 'EPICC-085-36', 'EPICC-085-45'))

agNames(map) <- c('B.1.1.7', 'B.1.351', 'B.1.429', 'B.1.617.2', 'P.1',
                  'C.37', 'B.1.621', 'B.1.526+E484K', 'D614G', 'BA.1',
                  'BA.1.1', 'BA.2', 'BA.2.12.1', 'BA.5', 'R.1')

agFill(map) <- c('#637939', '#e7ba52', '#9B9FD9', '#d18652', '#7b4173',
                 '#79c9c9', '#0096ad', '#9ab370', '#393b79', '#EF3737',
                 '#EF3737', '#d10fa2', '#d10fa2', '#F08DA5', '#f2f53b')

agSize(map) <- 18
srSize(map) <- 10
srOutlineWidth(map) <- 1
srFill(map) <- 'transparent'

srGroups(map) <- c('D614G convalescent', 'D614G convalescent', 'D614G convalescent', 'D614G convalescent',
                   'D614G convalescent', 'D614G convalescent', 'B.1.429 convalescent', 'B.1.429 convalescent',
                   'D614G convalescent', 'B.1.429 convalescent', 'B.1.1.7 convalescent', 'B.1.1.7 convalescent',
                   'B.1.1.7 convalescent', 'B.1.1.7 convalescent', 'B.1.1.7 convalescent', 'B.1.1.7 convalescent',
                   'C.37 convalescent', 'C.37 convalescent', 'C.37 convalescent', 'C.37 convalescent',
                   'B.1.1.7 convalescent', 'P.1 convalescent', 'C.37 convalescent', 'C.37 convalescent',
                   'C.37 convalescent', 'C.37 convalescent', 'C.37 convalescent', 'C.37 convalescent',
                   'P.1 convalescent', 'P.1 convalescent', 'B.1.351 convalescent', 'D614G convalescent',
                   'D614G convalescent', 'D614G convalescent', 'D614G convalescent', 'D614G convalescent',
                   'D614G convalescent', 'D614G convalescent', 'D614G convalescent', 'D614G convalescent',
                   'B.1.617.2 convalescent', 'B.1.617.2 convalescent', 'B.1.617.2 convalescent',
                   'B.1.617.2 convalescent', 'BA.1 convalescent', 'D614G convalescent', 'D614G convalescent',
                   'other', 'B.1.1.7 convalescent', 'P.1 convalescent', 'B.1.617.2 convalescent',
                   'B.1.617.2 convalescent', 'B.1.617.2 convalescent', 'B.1.617.2 convalescent',
                   'B.1.617.2 convalescent', 'BA.1 convalescent', 'B.1.429 convalescent', 'B.1.429 convalescent',
                   'B.1.617.2 convalescent', 'B.1.1.7 convalescent', 'B.1.617.2 convalescent', 'D614G convalescent',
                   'D614G convalescent', 'D614G convalescent', 'D614G convalescent', 'D614G convalescent',
                   'D614G convalescent', 'B.1.617.2 convalescent', 'B.1.617.2 convalescent', 'D614G convalescent',
                   'B.1.526+E484K convalescent', 'B.1.1.7 convalescent', 'B.1.617.2 convalescent', 'B.1.617.2 convalescent',
                   'B.1.351 convalescent', 'B.1.351 convalescent', 'B.1.351 convalescent', 'B.1.351 convalescent',
                   'B.1.351 convalescent', 'B.1.351 convalescent')

srOutline(map)[srGroups(map) == 'D614G convalescent'] <- '#333333'
srOutline(map)[srGroups(map) == 'B.1.429 convalescent'] <- '#9B9FD9'
srOutline(map)[srGroups(map) == 'B.1.1.7 convalescent'] <- '#637939'
srOutline(map)[srGroups(map) == 'C.37 convalescent'] <- '#79c9c9'
srOutline(map)[srGroups(map) == 'P.1 convalescent'] <- '#7b4173'
srOutline(map)[srGroups(map) == 'B.1.351 convalescent'] <- '#e7ba52'
srOutline(map)[srGroups(map) == 'B.1.617.2 convalescent'] <- '#d18652'
srOutline(map)[srGroups(map) == 'BA.1 convalescent'] <- '#EF3737'
srOutline(map)[srGroups(map) == 'B.1.526+E484K convalescent'] <- '#9ab370'
srOutline(map)[srGroups(map) == 'other'] <- 'white'

map <- optimizeMap(map, 2, 1000)

m <- read.acmap('data/individual_datasets/duke/maps/map_no_outliers.ace')

map <- realignMap(map, m)

map <- rotateMap(map, -34)

save.acmap(map, 'data/individual_datasets/fda/maps/map.ace')


