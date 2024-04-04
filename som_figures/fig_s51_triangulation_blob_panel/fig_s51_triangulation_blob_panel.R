
rm(list = ls())
library(Racmacs)

source('code/plotting/plot_maps.R')
source('code/data_generation/load_maps_for_comparison.R')

# Read in the maps
maps <- load_maps_for_comparison('', sr_group_standardisation = 'data/metadata/sr_group_standardisation_WT_conv_as_D614G_as614D.csv')

# Make maps with triangulation blobs with a stress limit of 1% of the map stress
sl = 0.01
gs = 0.1

duke <- triangulationBlobs(maps[[1]], stress_lim = mapStress(maps[[1]]) * sl, grid_spacing = gs)
maryland <- triangulationBlobs(maps[[2]], stress_lim = mapStress(maps[[2]]) * sl, grid_spacing = gs)
galveston <- triangulationBlobs(maps[[3]], stress_lim = mapStress(maps[[3]]) * sl, grid_spacing = gs)
emory <- triangulationBlobs(maps[[4]], stress_lim = mapStress(maps[[4]]) * sl, grid_spacing = gs)
madison_pooled <- triangulationBlobs(maps[[5]], stress_lim = mapStress(maps[[5]]) * sl, grid_spacing = gs)
madison_unpooled <- triangulationBlobs(maps[[6]], stress_lim = mapStress(maps[[6]]) * sl, grid_spacing = gs)
st_louis <- triangulationBlobs(maps[[7]], stress_lim = mapStress(maps[[7]]) * sl, grid_spacing = gs)
oxford <- triangulationBlobs(maps[[8]], stress_lim = mapStress(maps[[8]]) * sl, grid_spacing = gs)
mt_sinai_human <- triangulationBlobs(maps[[9]], stress_lim = mapStress(maps[[9]]) * sl, grid_spacing = gs)
emc <- triangulationBlobs(maps[[10]], stress_lim = mapStress(maps[[10]]) * sl, grid_spacing = gs)
innsbruck <- triangulationBlobs(maps[[11]], stress_lim = mapStress(maps[[11]]) * sl, grid_spacing = gs)
charite <- triangulationBlobs(maps[[12]], stress_lim = mapStress(maps[[12]]) * sl, grid_spacing = gs)
madison_frnt <- triangulationBlobs(maps[[13]], stress_lim = mapStress(maps[[13]]) * sl, grid_spacing = gs)
fda <- triangulationBlobs(maps[[14]], stress_lim = mapStress(maps[[14]]) * sl, grid_spacing = gs)
geneva <- triangulationBlobs(maps[[15]], stress_lim = mapStress(maps[[15]]) * sl, grid_spacing = gs)
amc <- triangulationBlobs(maps[[16]], stress_lim = mapStress(maps[[16]]) * sl, grid_spacing = gs)
emc_calu <- triangulationBlobs(maps[[17]], stress_lim = mapStress(maps[[17]]) * sl, grid_spacing = gs)
emc_vero <- triangulationBlobs(maps[[18]], stress_lim = mapStress(maps[[18]]) * sl, grid_spacing = gs)

# Do the plotting
png('som_figures/fig_s51_triangulation_blob_panel/fig_s51_triangulation_blob_panel.png', 25, 11, units = 'in', res=300, pointsize = 12)
layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 12, 13, 14, 15, 0, 0, 0, 0, 0, 0, 0, 0, 16, 17, 18, 0, 0, 0), ncol=6, byrow = T), heights=c(3, 3, 3, 0.5, 3))
par(mai=c(0.1, 0.1, 0.1, 0.1))
for (map in list(duke, fda, innsbruck, amc, emory, geneva, 
                 charite, madison_frnt, emc, emc_vero, emc_calu,
                 oxford, mt_sinai_human, st_louis, galveston,
                 maryland, madison_unpooled, madison_pooled)) {
  print(mapName(map))

  plot(map, xlim = c(-1, 9), ylim = c(-3, 3), cex = 1.3)

  text(9, 2.5, pretty_labels[mapName(map)], cex = 2, pos = 2, font = 2)
  text(9, -2.6, dataset_labels_sub[mapName(map)], cex = 1.8, pos = 2, font = 1)

}
dev.off()

