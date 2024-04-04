
rm(list = ls())

library(tidyverse)
library(Racmacs)

source('code/plotting/viewer_plotting.R')

# Read in the map
map <- read.acmap('data/merged_map/maps/merged_duplicated_antigens_only.ace')

# Run the dimensionality testing
mapDims <- dimensionTestMap(
  map,
  dimensions_to_test = 1:5,
  test_proportion = 0.1,
  minimum_column_basis = 'none',
  fixed_column_bases = colBases(map),
  number_of_optimizations = 1000,
  replicates_per_dimension = 100,
  options = list()
)

# Plot the figure
df <- data.frame(dimensions=c(1, 2, 3, 4, 5),
                 rmse=c(mapDims$mean_rmse_detectable))

ggplot(data=df, aes(x=dimensions, y=rmse)) +
  geom_line()+
  geom_point() +
  theme_bw() +
  xlab('Dimension') +
  ylab('Mean RMSE of detectable titers') +
  theme(strip.background = element_blank())
ggsave('som_figures/fig_s54_merged_map_dimensionality_test/fig_s54_merged_map_dimensionality_test_dimensionality_test.pdf',
       width = 4, height=3)

ggplot(data=df, aes(x=dimensions, y=rmse)) +
  geom_line()+
  geom_point() +
  theme_bw() +
  coord_cartesian(
    ylim = c(0, 2)
  ) +
  xlab('Dimension') +
  ylab('Mean RMSE of detectable titers') +
  theme(strip.background = element_blank())
ggsave('som_figures/fig_s54_merged_map_dimensionality_test/fig_s54_merged_map_dimensionality_test_dimensionality_test_axes.pdf',
       width = 4, height=3)

# Print the results
mapDims
# dimensions mean_rmse_detectable var_rmse_detectable mean_rmse_nondetectable var_rmse_nondetectable replicates
# 1          1             1.896587         0.009501841                2.048643            0.021550244        100
# 2          2             1.423079         0.009469993                1.606970            0.012100762        100
# 3          3             1.427599         0.007901530                1.497223            0.009468922        100
# 4          4             1.451486         0.013335105                1.454498            0.010434492        100
# 5          5             1.455702         0.013697523                1.441265            0.014199372        100


# Plot a 3D map
map3D <- read.acmap('data/merged_map/maps/merged_duplicated_antigens_only_3D.ace')

procrustes_2d_to_3d <- function(
    map2d, 
    map3d,
    link_col = 'black',
    point_opacity_2d = 1
) {
  
  map3d <- add_procrustes_grid(map3d, map2d)
  map3d <- add_2d_points(map3d, map2d, point_scaling = 0.15, point_opacity = point_opacity_2d)
  map3d <- link_2d_points(map3d, map2d, col = link_col)
  map3d
  
}

plot_map(
  map = procrustes_2d_to_3d(map, map3D),
  scaling_size = 0.15,
  sr_scaling_size = 0.1,
  alter_opacity = F,
  sera_opacity = 0.2,
  grid.col = NA
)

