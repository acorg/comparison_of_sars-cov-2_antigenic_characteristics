
rm(list = ls())

library(Racmacs)

source('code/metadata/common.R')

# Read in merged map for comparison
map <- read.acmap('data/merged_map/maps/merged_duplicated_antigens_only.ace')

srOutline(map) <- adjustcolor(srOutline(map), alpha.f = 0.5)

map <- rotateMap(map, 4)
map <- translateMap(map, c(1, 0.75))
srOutlineWidth(map) <- 0.5
srSize(map) <- 7

# Color variants in merged map by substitutions
agFill(map) <- ag_colors_info$ag_cols[match(agNames(map), ag_colors_info$ags)]

png('som_figures/fig_s60_merged_map_individual_map_cross_validation/fig_s60_merged_map_individual_map_cross_validation.png', 25, 14, units = 'in', res=300, pointsize = 12)
layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 12, 13, 14, 15, 0, 0, 0, 0, 0, 0, 0, 0, 16, 17, 18, 0, 0, 0), ncol=6, byrow = T), heights=c(3, 3, 3, 0.5, 3))
for (mapN in list('duke', 'fda', 'innsbruck', 'amc', 'emory', 'geneva', 
                  'charite', 'madison_frnt', 'emc_prnt', 'emc_vero', 'emc_calu',
                  'oxford', 'mt_sinai_human', 'st_louis', 'galveston',
                  'maryland', 'madison_unpooled', 'madison_pooled')) {
  
  print(mapN)
  m_ <- read.acmap(paste0('data/merged_map/maps/merged_duplicated_antigens_without_', mapN, '.ace'))
  m_ <- realignMap(m_, map)
  
  agFill(m_) <- ag_colors_info$ag_cols[match(agNames(m_), ag_colors_info$ags)]
  p <- procrustesMap(m_, map, sera = F)
  
  par(mar = c(0.3, 0, 0.3, 0))
  plot(p, xlim = c(1, 9), ylim = c(-2, 4), fill.alpha = 0.9, plot_labels = FALSE,
       outline.alpha = 1, grid.col = '#cfcfcf', plot_stress = FALSE,
       grid.margin.col = 'black', show_error_lines = FALSE, cex = 1.5)
  
  text(9, 3.5, pretty_labels[mapN], cex = 2, pos = 2, font = 2)
}

dev.off()
