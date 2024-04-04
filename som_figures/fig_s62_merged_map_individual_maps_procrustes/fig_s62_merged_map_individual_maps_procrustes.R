
rm(list = ls())

library(Racmacs)

source('code/data_generation/load_maps_for_comparison.R')
source('code/metadata/common.R')

# Get individual maps for comparison
maps <- load_maps_for_comparison('', sr_group_standardisation = 'data/metadata/sr_group_standardisation_WT_conv_as_D614G_as614D.csv')

# Read in merged map for comparison
map <- read.acmap('data/merged_map/maps/merged_duplicated_antigens_only.ace')

# Change antigen colors
for(m_ in names(maps)) {
  agFill(maps[[m_]]) <- ag_colors_info$ag_cols[match(agNames(maps[[m_]]), ag_colors_info$ags, nomatch = 29)]
}

agFill(map) <- ag_colors_info$ag_cols[match(agNames(map), ag_colors_info$ags)]

# Style the merged map
srOutline(map) <- adjustcolor(srOutline(map), alpha.f = 0.5)

map <- rotateMap(map, 4)
map <- translateMap(map, c(1, 0.75))
srOutlineWidth(map) <- 0.5
srSize(map) <- 7


# Figure colored by substitution

png('som_figures/fig_s62_merged_map_individual_maps_procrustes/fig_s62_merged_map_individual_maps_procrustes.png', 30, 15, units = 'in', res=300, pointsize = 12)
layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 12, 13, 14, 15, 0, 0, 0, 0, 0, 0, 0, 0, 16, 17, 18, 0, 0, 0), ncol=6, byrow = T), heights=c(3, 3, 3, 0.5, 3))
for (mapName in list('duke', 'fda', 'innsbruck', 'amc', 'emory', 'geneva',
                     'charite', 'madison_frnt', 'emc_prnt', 'emc_vero', 'emc_calu',
                     'oxford', 'mt_sinai_human', 'st_louis', 'galveston',
                     'maryland', 'madison_unpooled', 'madison_pooled')) {

  m_ <- maps[[mapName]]
  
  p <- procrustesMap(map, m_, sera = F)

  par(mar = c(0.3, 0, 0.3, 0))
  plot(p, xlim = c(1, 10), ylim = c(-2, 4), fill.alpha = 1, plot_labels = FALSE,
       outline.alpha = 1, grid.col = '#cfcfcf', plot_stress = FALSE,
       grid.margin.col = 'black', show_error_lines = FALSE, cex = 1.5)

  text(10, 3.5, prett_labels[mapName(m_)], cex = 2, pos = 2, font = 2)
}

dev.off()
