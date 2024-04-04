
rm(list = ls())
library(Racmacs)
set.seed(100)

map <- read.acmap('data/merged_map/maps/merged_duplicated_antigens_only.ace')

# Resampling antigens and sera
mapBootAS <- bootstrapMap(
  map,
  'resample',
  bootstrap_repeats = 100,
  bootstrap_ags = TRUE,
  bootstrap_sr = TRUE,
  reoptimize = TRUE,
  optimizations_per_repeat = 1000,
  ag_noise_sd = 0.4,
  titer_noise_sd = 0.62,
  options = list()
)

mapBootASBlobs <- bootstrapBlobs(mapBootAS, conf.level = 0.68, smoothing = 2, gridspacing = 0.05)

# Set plotting limits
plotlims <- Racmacs:::mapPlotLims(mapBootASBlobs, sera = F)

# Plot the figure
png('som_figures/fig_s57_merged_map_resample_bootstrap/fig_s57_merged_map_resample_bootstrap.png', width = 6, height = 4.7, units = 'in', res=300, pointsize = 18)
plot(mapBootASBlobs, xlim = plotlims$xlim, ylim = plotlims$ylim, fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
     grid.col = '#cfcfcf', grid.margin.col='#7d7d7d', cex=0.7)
dev.off()

