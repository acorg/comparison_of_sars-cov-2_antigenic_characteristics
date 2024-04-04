
rm(list = ls())

library(Racmacs)

source('code/metadata/common.R')

# Read in merged map for comparison
map <- read.acmap('data/merged_map/maps/merged_duplicated_antigens_only.ace')

# Style the merged map
srOutline(map) <- adjustcolor(srOutline(map), alpha.f = 0.5)
map <- rotateMap(map, 4)
map <- translateMap(map, c(1, 0.75))
srOutlineWidth(map) <- 0.5
srSize(map) <- 7

# Figure colored by substitution
agFill(map) <- ag_colors_info$ag_cols[match(agNames(map), ag_colors_info$ags)]


# Plot procrustes for map without cpe map
png('som_figures/fig_s61_merged_map_cpe_maps_cross_validation/fig_s61_merged_map_cpe_maps_cross_validation.png', 12, 5, units = 'in', res=300, pointsize = 12)
layout(matrix(c(1, 2), ncol=2, byrow = T))

map <- translateMap(map, c(0.2, 0.2))

m_ <- read.acmap(paste0('data/merged_map/maps/merged_duplicated_antigens_without_cpe.ace'))
m_ <- realignMap(m_, map)
agFill(m_) <- ag_colors_info$ag_cols[match(agNames(m_), ag_colors_info$ags)]

p_cpe <- procrustesMap(map, m_, sera = F)
p_m <- procrustesMap(m_, map, sera = F)

par(mar = c(0.3, 0, 0.3, 0))
plot(p_cpe, xlim = c(1, 9), ylim = c(-2, 4), fill.alpha = 0.9, plot_labels = FALSE,
     outline.alpha = 1, grid.col = '#cfcfcf', plot_stress = FALSE,
     grid.margin.col = 'black', show_error_lines = FALSE, cex = 1)

text(1.5, 3.5, 'A', cex = 2)

par(mar = c(0.3, 0, 0.3, 0))
plot(p_m, xlim = c(1, 9), ylim = c(-2, 4), fill.alpha = 0.9, plot_labels = FALSE,
     outline.alpha = 1, grid.col = '#cfcfcf', plot_stress = FALSE,
     grid.margin.col = 'black', show_error_lines = FALSE, cex = 1)

text(1.5, 3.5, 'B', cex = 2)

dev.off()

