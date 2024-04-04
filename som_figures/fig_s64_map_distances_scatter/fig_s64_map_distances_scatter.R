
rm(list = ls())

library(Racmacs)
library(tidyverse)
library(patchwork)

source('code/data_generation/load_maps_for_comparison.R')
source('code/metadata/common.R')

# Read in the individual maps
maps <- load_maps_for_comparison('', sr_group_standardisation = 'data/metadata/sr_group_standardisation_WT_conv_as_D614G_as614D.csv')

# Read in merged map for comparison
map <- read.acmap('data/merged_map/maps/merged_duplicated_antigens_only.ace')

# Get the map distances from the merged map
mapDists <- mapDistances(map)
colnames(mapDists) <- srNames(map)
rownames(mapDists) <- agNames(map)


# Plot map vs map distances
toplot <- lapply(list('duke', 'fda', 'innsbruck', 'amc', 'emory', 'geneva',
                      'charite', 'madison_frnt', 'emc_prnt', 'emc_vero', 'emc_calu',
                      'oxford', 'mt_sinai_human', 'st_louis', 'galveston',
                      'maryland', 'madison_unpooled', 'madison_pooled'), function(mapName) {
                        
                        map <- maps[[mapName]]
                        print(mapName(map))
                        
                        map <- removeAntigens(map, agNames(map)[!(agNames(map) %in% duplicated_ags)])
                        tDists <- mapDistances(map)
                        colnames(tDists) <- srNames(map)
                        rownames(tDists) <- agNames(map)
                        
                        mapDistsMod <- mapDists[rownames(mapDists) %in% agNames(map), colnames(mapDists) %in% srNames(map)]
                        
                        ggplot(
                          data = tibble(
                            x = as.numeric(c(tDists)), 
                            y = as.numeric(c(mapDistsMod)), 
                            Variants = c(matrix(unlist(rep(rownames(tDists), length(colnames(tDists)))), ncol = length(colnames(tDists)), nrow = length(rownames(tDists))))
                            ),
                          aes(
                            x = x,
                            y = y,
                            color = Variants
                          )
                        ) +
                          geom_abline(
                            slope = 1,
                            intercept = 0,
                            color = 'grey',
                            linetype = 'dashed'
                          ) +
                          geom_point() +
                          scale_color_manual(
                            values = ag_colors
                          ) +
                          theme_bw() +
                          coord_cartesian(
                            xlim = c(0, 11),
                            ylim = c(0, 11)
                          ) +
                          labs(
                            title = pretty_labels[mapName(map)],
                            x = paste0('map distances (', pretty_labels[mapName(map)], ')'),
                            y = 'map distances (merged map)'
                          ) +
                          guides(color = 'none')
                      })

gp <- wrap_plots(toplot, ncol = 3)
ggsave('som_figures/fig_s64_map_distances_scatter/fig_s64_map_distances_scatter.png', plot = gp, width = 12, height = 24)
ggsave('som_figures/fig_s64_map_distances_scatter/fig_s64_map_distances_scatter.pdf', plot = gp, width = 12, height = 24)

