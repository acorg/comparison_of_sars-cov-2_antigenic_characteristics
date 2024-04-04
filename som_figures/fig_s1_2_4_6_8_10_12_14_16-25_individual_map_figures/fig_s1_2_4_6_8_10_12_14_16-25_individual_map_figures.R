
rm(list = ls())

library(Racmacs)
library(tidyverse)

source('code/data_generation/load_maps_for_comparison.R')
source('code/functions/diagnostics.R')

test_dimensions <- function(map) {
  # Run the dimensionality testing
  mapDims <- dimensionTestMap(
    map,
    dimensions_to_test = 1:4,
    test_proportion = 0.1,
    minimum_column_basis = 'none',
    fixed_column_bases = colBases(map),
    number_of_optimizations = 500,
    replicates_per_dimension = 500,
    options = list()
  )
  
  # Plot the figure
  df <- data.frame(dimensions=c(1, 2, 3, 4),
                   rmse=c(mapDims$mean_rmse_detectable))
  
  ggplot(data=df, aes(x=dimensions, y=rmse)) +
    geom_line()+
    geom_point() +
    theme_bw() +
    xlab('Dimension') +
    ylab('Mean RMSE of detectable titers') +
    theme(strip.background = element_blank()) -> gp

  gp
}

map_vs_table <- function(map, filename) {
  # plot map versus table distances
  residual_data <- residualErrorTable(map)
  ggplot(
    data = residual_data) +
    geom_point(
      mapping = aes(
        y     = measured_logtiter_upper,
        x     = predicted_logtiter,
      ),
      alpha = 0.4
    ) +
    theme_bw() +
    xlim(1, 10) + 
    ylim(1, 10) +
    ylab('Measured log titers') + 
    xlab('Fitted log titers') + 
    geom_abline(intercept = 0, lty=2, col = 'grey') + 
    theme(text = element_text(size=15)) -> gp
  
  ggsave(filename, width=5.5, height=5)
}
  

# Read in data
maps <- load_maps_for_comparison('', sr_group_standardisation = 'data/metadata/sr_group_standardisation_WT_conv_as_D614G_as614D.csv')

duke <- maps[[1]]
emory <- maps[[4]]
fda <- maps[[14]]
innsbruck <- maps[[11]]
oxford <- maps[[8]]
mt_sinai_human <- maps[[9]]
amc <- maps[[16]]
geneva <- maps[[15]]
madison_pooled <- maps[[5]]
charite <- maps[[12]]
emc <- maps[[10]]
madison_unpooled <- maps[[6]]
emc_vero <- maps[[18]]
emc_calu <- maps[[17]]
galveston <- maps[[3]]
madison_frnt <- maps[[13]]
maryland <- maps[[2]]
wustl <- maps[[7]]


## Make individual map figures

# Plot duke map
ggplot(duke) -> gp
gp +
  ggrepel::geom_label_repel(
    data = filter(gp$data, type == 'AG'),
    mapping = aes(
      x = x,
      y = y,
      label = text,
      point.size = size
    ),
    size = 4,
    min.segment.length = 0,
    label.padding = 0.2,
    label.size = NA,
    fill = NA
  )
ggsave('som_figures/fig_s1_2_4_6_8_10_12_14_16-25_individual_map_figures/fig_s1_duke_map.png', width = 9, height = 9.5)

# Plot emory
ggplot(emory) -> gp
gp +
  ggrepel::geom_label_repel(
    data = filter(gp$data, type == 'AG'),
    mapping = aes(
      x = x,
      y = y,
      label = text,
      point.size = size
    ),
    size = 4,
    min.segment.length = 0,
    label.padding = 0.2,
    label.size = NA,
    fill = NA
  )
ggsave('som_figures/fig_s1_2_4_6_8_10_12_14_16-25_individual_map_figures/fig_s2_emory_map.png', width = 8, height = 7)

# Plot FDA
ggplot(fda) -> gp
gp +
  ggrepel::geom_label_repel(
    data = filter(gp$data, type == 'AG'),
    mapping = aes(
      x = x,
      y = y,
      label = text,
      point.size = size
    ),
    size = 4,
    min.segment.length = 0,
    label.padding = 0.2,
    label.size = NA,
    fill = NA
  )
ggsave('som_figures/fig_s1_2_4_6_8_10_12_14_16-25_individual_map_figures/fig_s4_fda_map.png', width = 6, height = 5)

# Plot Innsbruck
ggplot(innsbruck) -> gp
gp +
  ggrepel::geom_label_repel(
    data = filter(gp$data, type == 'AG'),
    mapping = aes(
      x = x,
      y = y,
      label = text,
      point.size = size
    ),
    size = 4,
    min.segment.length = 0,
    label.padding = 0.2,
    label.size = NA,
    fill = NA
  )
ggsave('som_figures/fig_s1_2_4_6_8_10_12_14_16-25_individual_map_figures/fig_s6_innsbruck_map.png', width = 6, height = 6.5)

# Plot Oxford
ggplot(oxford) -> gp
gp +
  ggrepel::geom_label_repel(
    data = filter(gp$data, type == 'AG'),
    mapping = aes(
      x = x,
      y = y,
      label = text,
      point.size = size
    ),
    size = 4,
    min.segment.length = 0,
    label.padding = 0.2,
    label.size = NA,
    fill = NA
  )
ggsave('som_figures/fig_s1_2_4_6_8_10_12_14_16-25_individual_map_figures/fig_s8_oxford_map.png', width = 6, height = 6.5)

# Plot Mt. Sinai
ggplot(mt_sinai_human) -> gp
gp +
  ggrepel::geom_label_repel(
    data = filter(gp$data, type == 'AG'),
    mapping = aes(
      x = x,
      y = y,
      label = text,
      point.size = size
    ),
    size = 4,
    min.segment.length = 0,
    label.padding = 0.2,
    label.size = NA,
    fill = NA
  )
ggsave('som_figures/fig_s1_2_4_6_8_10_12_14_16-25_individual_map_figures/fig_s10_mt_sinai_map.png', width = 4.5, height = 5)


# Plot AMC
ggplot(amc) -> gp
gp +
  ggrepel::geom_label_repel(
    data = filter(gp$data, type == 'AG'),
    mapping = aes(
      x = x,
      y = y,
      label = text,
      point.size = size
    ),
    size = 4,
    min.segment.length = 0,
    label.padding = 0.2,
    label.size = NA,
    fill = NA
  )
ggsave('som_figures/fig_s1_2_4_6_8_10_12_14_16-25_individual_map_figures/fig_s12_amc_map.png', width = 4.5, height = 4)

# Plot Geneva
ggplot(geneva) -> gp
gp +
  ggrepel::geom_label_repel(
    data = filter(gp$data, type == 'AG'),
    mapping = aes(
      x = x,
      y = y,
      label = text,
      point.size = size
    ),
    size = 4,
    min.segment.length = 0,
    label.padding = 0.2,
    label.size = NA,
    fill = NA
  )
ggsave('som_figures/fig_s1_2_4_6_8_10_12_14_16-25_individual_map_figures/fig_s14_geneva_map.png', width = 5, height = 4.5)

# Plot Madison (pooled)
ggplot(madison_pooled) -> gp
gp +
  ggrepel::geom_label_repel(
    data = filter(gp$data, type == 'AG'),
    mapping = aes(
      x = x,
      y = y,
      label = text,
      point.size = size
    ),
    size = 4,
    min.segment.length = 0,
    label.padding = 0.2,
    label.size = NA,
    fill = NA
  )
ggsave('som_figures/fig_s1_2_4_6_8_10_12_14_16-25_individual_map_figures/fig_s16_madison_pooled_map.png', width = 5, height = 5.5)

# Plot Charite
ggplot(charite) -> gp
gp +
  ggrepel::geom_label_repel(
    data = filter(gp$data, type == 'AG'),
    mapping = aes(
      x = x,
      y = y,
      label = text,
      point.size = size
    ),
    size = 4,
    min.segment.length = 0,
    label.padding = 0.2,
    label.size = NA,
    fill = NA
  )
ggsave('som_figures/fig_s1_2_4_6_8_10_12_14_16-25_individual_map_figures/fig_s17_charite_map.png', width = 5, height = 4.5)

# Plot EMC PRNT
ggplot(emc) -> gp
gp +
  ggrepel::geom_label_repel(
    data = filter(gp$data, type == 'AG'),
    mapping = aes(
      x = x,
      y = y,
      label = text,
      point.size = size
    ),
    size = 4,
    min.segment.length = 0,
    label.padding = 0.2,
    label.size = NA,
    fill = NA
  )
ggsave('som_figures/fig_s1_2_4_6_8_10_12_14_16-25_individual_map_figures/fig_s18_emc_map.png', width = 7, height = 5)

# Plot Madison (unpooled)
ggplot(madison_unpooled) -> gp
gp +
  ggrepel::geom_label_repel(
    data = filter(gp$data, type == 'AG'),
    mapping = aes(
      x = x,
      y = y,
      label = text,
      point.size = size
    ),
    size = 4,
    min.segment.length = 0,
    label.padding = 0.2,
    label.size = NA,
    fill = NA
  )
ggsave('som_figures/fig_s1_2_4_6_8_10_12_14_16-25_individual_map_figures/fig_s19_madison_unpooled_map.png', width = 4.5, height = 5)

# Plot EMC VeroE6
ggplot(emc_vero) -> gp
gp +
  ggrepel::geom_label_repel(
    data = filter(gp$data, type == 'AG'),
    mapping = aes(
      x = x,
      y = y,
      label = text,
      point.size = size
    ),
    size = 4,
    min.segment.length = 0,
    label.padding = 0.2,
    label.size = NA,
    fill = NA
  )
ggsave('som_figures/fig_s1_2_4_6_8_10_12_14_16-25_individual_map_figures/fig_s20_emc_vero_map.png', width = 6, height = 4)

# Plot EMC Calu-3
ggplot(emc_calu) -> gp
gp +
  ggrepel::geom_label_repel(
    data = filter(gp$data, type == 'AG'),
    mapping = aes(
      x = x,
      y = y,
      label = text,
      point.size = size
    ),
    size = 4,
    min.segment.length = 0,
    label.padding = 0.2,
    label.size = NA,
    fill = NA
  )
ggsave('som_figures/fig_s1_2_4_6_8_10_12_14_16-25_individual_map_figures/fig_s21_emc_calu_map.png', width = 5.5, height = 4)

# Plot Galveston
ggplot(galveston) -> gp
gp +
  ggrepel::geom_label_repel(
    data = filter(gp$data, type == 'AG'),
    mapping = aes(
      x = x,
      y = y,
      label = text,
      point.size = size
    ),
    size = 4,
    min.segment.length = 0,
    label.padding = 0.2,
    label.size = NA,
    fill = NA
  )
ggsave('som_figures/fig_s1_2_4_6_8_10_12_14_16-25_individual_map_figures/fig_s22_galveston_map.png', width = 4, height = 4.5)

# Plot Madison (FRNT)
ggplot(madison_frnt) -> gp
gp +
  ggrepel::geom_label_repel(
    data = filter(gp$data, type == 'AG'),
    mapping = aes(
      x = x,
      y = y,
      label = text,
      point.size = size
    ),
    size = 4,
    min.segment.length = 0,
    label.padding = 0.2,
    label.size = NA,
    fill = NA
  )
ggsave('som_figures/fig_s1_2_4_6_8_10_12_14_16-25_individual_map_figures/fig_s23_madison_frnt_map.png', width = 6, height = 4.5)

# Plot Maryland
ggplot(maryland) -> gp
gp +
  ggrepel::geom_label_repel(
    data = filter(gp$data, type == 'AG'),
    mapping = aes(
      x = x,
      y = y,
      label = text,
      point.size = size
    ),
    size = 4,
    min.segment.length = 0,
    label.padding = 0.2,
    label.size = NA,
    fill = NA
  )
ggsave('som_figures/fig_s1_2_4_6_8_10_12_14_16-25_individual_map_figures/fig_s24_maryland_map.png', width = 6, height = 5)

# Plot WUSTL
ggplot(wustl) -> gp
gp +
  ggrepel::geom_label_repel(
    data = filter(gp$data, type == 'AG'),
    mapping = aes(
      x = x,
      y = y,
      label = text,
      point.size = size
    ),
    size = 4,
    min.segment.length = 0,
    label.padding = 0.2,
    label.size = NA,
    fill = NA
  )
ggsave('som_figures/fig_s1_2_4_6_8_10_12_14_16-25_individual_map_figures/fig_s25_wustl_map.png', width = 3.5, height = 5)

## Plot map vs table distances
map_vs_table(emory, 'som_figures/fig_s1_2_4_6_8_10_12_14_16-25_individual_map_figures/fig_s2_emory_map_vs_table.png')
map_vs_table(oxford, 'som_figures/fig_s1_2_4_6_8_10_12_14_16-25_individual_map_figures/fig_s8_oxford_map_vs_table.png')
map_vs_table(mt_sinai_human, 'som_figures/fig_s1_2_4_6_8_10_12_14_16-25_individual_map_figures/fig_s10_mt_sinai_map_vs_table.png')
map_vs_table(madison_pooled, 'som_figures/fig_s1_2_4_6_8_10_12_14_16-25_individual_map_figures/fig_s16_madison_pooled_map_vs_table.png')
map_vs_table(madison_unpooled, 'som_figures/fig_s1_2_4_6_8_10_12_14_16-25_individual_map_figures/fig_s19_madison_unpooled_map_vs_table.png')
map_vs_table(galveston, 'som_figures/fig_s1_2_4_6_8_10_12_14_16-25_individual_map_figures/fig_s22_galveston_map_vs_table.png')
map_vs_table(madison_frnt, 'som_figures/fig_s1_2_4_6_8_10_12_14_16-25_individual_map_figures/fig_s23_madison_frnt_map_vs_table.png')
map_vs_table(maryland, 'som_figures/fig_s1_2_4_6_8_10_12_14_16-25_individual_map_figures/fig_s24_maryland_map_vs_table.png')
map_vs_table(wustl, 'som_figures/fig_s1_2_4_6_8_10_12_14_16-25_individual_map_figures/fig_s25_wustl_map_vs_table.png')

## Make 3D maps
emory3D <- optimizeMap(emory, 3, 500)
p <- procrustesMap(emory3D, emory, sera = FALSE)
view(p)

oxford3D <- optimizeMap(oxford, 3, 500)
p <- procrustesMap(oxford3D, oxford, sera = FALSE)
view(p)

mt_sinai_human3D <- optimizeMap(mt_sinai_human, 3, 500)
p <- procrustesMap(mt_sinai_human3D, mt_sinai_human, sera = FALSE)
view(p)

madison_pooled3D <- optimizeMap(madison_pooled, 3, 500)
p <- procrustesMap(madison_pooled3D, madison_pooled, sera = FALSE)
view(p)

madison_unpooled3D <- optimizeMap(madison_unpooled, 3, 500)
p <- procrustesMap(madison_unpooled3D, madison_unpooled, sera = FALSE)
view(p)

galveston3D <- optimizeMap(galveston, 3, 500)
p <- procrustesMap(galveston3D, galveston, sera = FALSE)
view(p)

madison_frnt3D <- optimizeMap(madison_frnt, 3, 500)
p <- procrustesMap(madison_frnt3D, madison_frnt, sera = FALSE)
view(p)

maryland3D <- optimizeMap(maryland, 3, 500)
p <- procrustesMap(maryland3D, maryland, sera = FALSE)
view(p)

wustl3D <- optimizeMap(wustl, 3, 500)
p <- procrustesMap(wustl3D, wustl, sera = FALSE)
view(p)

## Do dimensionality tests
gp_emory <- test_dimensions(emory)
ggsave('som_figures/fig_s1_2_4_6_8_10_12_14_16-25_individual_map_figures/fig_s2_emory_dimensionality_test.png', gp_emory, width = 4, height=3)

gp_oxford <- test_dimensions(oxford)
ggsave('som_figures/fig_s1_2_4_6_8_10_12_14_16-25_individual_map_figures/fig_s8_oxford_dimensionality_test.png', gp_oxford, width = 4, height=3)

gp_mt_sinai_human <- test_dimensions(mt_sinai_human)
ggsave('som_figures/fig_s1_2_4_6_8_10_12_14_16-25_individual_map_figures/fig_s10_mt_sinai_dimensionality_test.png', gp_mt_sinai_human, width = 4, height=3)

gp_madison_pooled <- test_dimensions(madison_pooled)
ggsave('som_figures/fig_s1_2_4_6_8_10_12_14_16-25_individual_map_figures/fig_s16_madison_pooled_dimensionality_test.png', gp_madison_pooled, width = 4, height=3)

gp_madison_unpooled <- test_dimensions(madison_unpooled)
ggsave('som_figures/fig_s1_2_4_6_8_10_12_14_16-25_individual_map_figures/fig_s19_madison_unpooled_dimensionality_test.png', gp_madison_unpooled, width = 4, height=3)

gp_galveston <- test_dimensions(galveston)
ggsave('som_figures/fig_s1_2_4_6_8_10_12_14_16-25_individual_map_figures/fig_s22_galveston_dimensionality_test.png', gp_galveston, width = 4, height=3)

gp_madison_frnt <- test_dimensions(madison_frnt)
ggsave('som_figures/fig_s1_2_4_6_8_10_12_14_16-25_individual_map_figures/fig_s23_madison_frnt_dimensionality_test.png', gp_madison_frnt, width = 4, height=3)

gp_maryland <- test_dimensions(maryland)
ggsave('som_figures/fig_s1_2_4_6_8_10_12_14_16-25_individual_map_figures/fig_s24_maryland_dimensionality_test.png', gp_maryland, width = 4, height=3)

gp_wustl <- test_dimensions(wustl)
ggsave('som_figures/fig_s1_2_4_6_8_10_12_14_16-25_individual_map_figures/fig_s25_wustl_dimensionality_test.png', gp_wustl, width = 4, height=3)

