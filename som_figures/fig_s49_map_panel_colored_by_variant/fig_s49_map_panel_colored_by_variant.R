
rm(list = ls())

library(Racmacs)

source('code/plotting/plot_maps.R')
source('code/data_generation/load_maps_for_comparison.R')

maps <- load_maps_for_comparison('', sr_group_standardisation = 'data/metadata/sr_group_standardisation_WT_conv_as_D614G_as614D.csv')

merged <- read.acmap('data/merged_map/maps/merged-duplicated-antigens-only.ace')

innsbruck_ <- translateMap(maps[[11]], c(0.4, -1))
duke_ <- translateMap(maps[[1]], c(-0.5, -0.5))
emory_ <- translateMap(maps[[4]], c(0, -0.9))
geneva_ <- translateMap(maps[[15]], c(0, -0.7))
madison_frnt_ <- translateMap(maps[[13]], c(0, -0.8))
emc_vero_ <- translateMap(maps[[18]], c(-0.2, -0.75))
emc_calu_ <- translateMap(maps[[17]], c(0, -0.6))
oxford_ <- translateMap(maps[[8]], c(0.2, -0.25))
mt_sinai_human_ <- translateMap(maps[[9]], c(0.1, -0.7))
st_louis_ <- translateMap(maps[[7]], c(0, -0.5))
madison_unpooled_ <- translateMap(maps[[6]], c(-0.3, -0.6))
charite_ <- translateMap(maps[[12]], c(-0.1, 0.6))
emc_ <- translateMap(maps[[10]], c(0, -0.2))
fda_ <- rotateMap(maps[[14]], 0.5)
fda_ <- translateMap(fda_, c(0.3, 0))
maryland_ <- translateMap(maps[[2]], c(0.1, 0))
galveston_ <- translateMap(maps[[3]], c(0.1, 0))
madison_pooled_ <- maps[[5]]
amc_ <- maps[[16]]

png('som_figures/fig_s49_map_panel_colored_by_variant/fig_s49_map_panel_colored_by_variant.png', 23, 10, units = 'in', res=300, pointsize = 12)
layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 12, 13, 14, 15, 0, 0, 0, 0, 0, 0, 0, 0, 16, 17, 18, 0, 0, 0), ncol=6, byrow = T), heights=c(3, 3, 3, 0.5, 3))
for (map_ in list(duke_, fda_, innsbruck_, amc_, emory_, geneva_, 
                  charite_, madison_frnt_, emc_prnt_, emc_vero_, emc_calu_,
                  oxford_, mt_sinai_human_, st_louis_, galveston_,
                  maryland_, madison_unpooled_, madison_pooled_)) {
  
  print(mapName(map_))
  
  map <- procrustesMap(map_, merged, sera = F)
  
  wt_ag <- unname(wt_ags[mapName(map)])
  agNames(map)[agNames(map) == 'BA.4/BA.5'] <- 'BA.5'
  specials <- c(wt_ag, 'B.1.351', 'B.1.617.2', 'BA.1', 'BA.2', 'BA.5')
  
  to_grey <- agNames(map)[! agNames(map) %in% specials]
  
  greyed_map(map, to_grey, specials, 0.5, xlim = c(-1, 9), ylim = c(-3, 3), cex = 0.8, procrustes = T)  #cex = 1.3, 1
  
  coords <- agCoords(map)
  
  lty = 1
  a = 0.3
  
  lines(c(coords[wt_ag, 1], coords['B.1.617.2', 1]), c(coords[wt_ag, 2], coords['B.1.617.2', 2]), col = adjustcolor('black', alpha=a), lty = lty, type = 'l', lwd = 4)
  if ('B.1.351' %in% agNames(map)) {
    lines(c(coords[wt_ag, 1], coords['B.1.351', 1]), c(coords[wt_ag, 2], coords['B.1.351', 2]), col = adjustcolor('black', alpha=a), lty = lty, type = 'l', lwd = 4)
  }
  if ('BA.1' %in% agNames(map) & ('BA.2' %in% agNames(map))){
    lines(c(coords['BA.2', 1], coords['BA.1', 1]), c(coords['BA.2', 2], coords['BA.1', 2]), col = adjustcolor('black', alpha=a), lty = lty, type = 'l', lwd = 4)
  }
  if ('BA.1' %in% agNames(map) & ('BA.5' %in% agNames(map))){
    lines(c(coords['BA.5', 1], coords['BA.1', 1]), c(coords['BA.5', 2], coords['BA.1', 2]), col = adjustcolor('black', alpha=a), lty = lty, type = 'l', lwd = 4)
  }
  if ('BA.1' %in% agNames(map)) {
    lines(c(coords[wt_ag, 1], coords['BA.1', 1]), c(coords[wt_ag, 2], coords['BA.1', 2]), col = adjustcolor('black', alpha=a), lty = lty, type = 'l', lwd = 4)
  }
  
  text(9, 2.45, pretty_labels[mapName(map)], cex = 2.2, pos = 2, font = 2)
  text(9, -2.6, dataset_labels_sub[mapName(map)], cex = 2.1, pos = 2, font = 1)
  
}
dev.off()


