
rm(list = ls())

library(tidyverse)
library(titertools)

source('code/metadata/common.R')
source('code/plotting/map_plotstyles.R')
source('code/plotting/plot_titers.R')
source('code/functions/gmts.R')

# Read in data
map_info <- readRDS('data/titer_analyses/titer_magnitude/dataset_magnitude_effect/dataset_magnitude_effect_map_info_subset.rds')

levels(map_info$sr_group)[10] <- 'mRNA-1273 or BNT162b2 vaccinated'


# Serum groups and antigens to consider
sr_groups <- c('mRNA-1273 or BNT162b2 vaccinated', 'D614G convalescent', 'B.1.1.7 convalescent', 
               'B.1.351 convalescent', 'P.1 convalescent', 'B.1.617.2 convalescent',
               'BA.1 convalescent')


map_info %>%
  filter(sr_group %in% sr_groups) %>%
  filter(ag_name %in% duplicated_ags) %>%
  group_by(
    map, sr_group, ag_name
  ) %>%
  summarise_gmts() %>%
  filter(!is.na(gmt)) -> map_gmts

# Calculate the gmts for all serum groups/antigens that should be plotted across all maps
map_gmts %>%
  group_by(
    sr_group, ag_name
  ) %>%
  summarise(gmt=median(gmt), gmt_upper=max(gmt), gmt_lower=min(gmt)) %>%
  ungroup() -> sr_group_gmts

# Order the serum groups, add information about animal model
map_gmts %>%
  mutate(
    sr_group = factor(sr_group, levels = sr_groups),
    map = factor(map, levels = animal_dataset_order),
    animal =  factor(unlist(lapply(as.character(map), function(i) model_animal[i,]), use.names=TRUE), levels = c('human', 'hamster', 'mouse'))
  ) -> map_gmts

# Rename columns for plotting
names(sr_group_gmts)[names(sr_group_gmts) == 'gmt'] <- 'estimate'
names(sr_group_gmts)[names(sr_group_gmts) == 'gmt_lower'] <- 'estimate_lower'
names(sr_group_gmts)[names(sr_group_gmts) == 'gmt_upper'] <- 'estimate_upper'

# Make the plot coloured by map
plot_adjusted_titer_comparison_gmt(map_gmts, sr_group_gmts, colors = map_colors, shapes = map_shapes,
                                   included_maps = unique(map_info$map), bar_alpha = 1, ag_order = duplicated_ags, 
                                   highlighted_map = c()) -> gp
gp <- gp + theme(
  panel.background = element_rect(fill = 'white', color = 'white'),
  strip.text = element_text(size=25),
  axis.text.x = element_text(size=25),
  axis.text.y = element_text(size=20),
  legend.text = element_text(size=20)
  ) +
  scale_y_titer(threshold = '<5') +
  scale_y_continuous(
    name = 'Titer',
    minor_breaks = NULL,
    breaks = seq(-2, 14),
    labels = c('<5', '', '10', '', '40', '', '160', '', '640', '', '2560', '', '10240', '', '40960', '', '163840')
  ) +
  coord_cartesian(ylim = c(-3, 14)) +
  scale_color_manual(
    values = map_colors,
    labels = unlist(lapply(animal_dataset_order, function(i) pretty_labels[i]),
                    use.names=TRUE)
  ) +
  labs(
    x = ''
  )
ggsave('som_figures/fig_s26-28_titer_magnitude/fig_s26_titer_magnitude_by_dataset.png', plot = gp, width = 22, height = 25)
ggsave('som_figures/fig_s26-28_titer_magnitude/fig_s26_titer_magnitude_by_dataset.pdf', plot = gp, width = 22, height = 25)


# Make the plot colored by animal model
plot_adjusted_titer_comparison_gmt(map_gmts, sr_group_gmts, colors = animal_colors, shapes = map_shapes,
                                   included_maps = unique(map_info$map), bar_alpha = 1, ag_order = duplicated_ags, 
                                   highlighted_map = c()) -> gp
gp <- gp + theme(
  panel.background = element_rect(fill = 'white', color = 'white'),
  strip.text = element_text(size=25),
  axis.text.x = element_text(size=25),
  axis.text.y = element_text(size=20),
  legend.text=element_text(size=20)
) +
  scale_y_titer(threshold = '<5') +
  scale_y_continuous(
    name = 'Titer',
    minor_breaks = NULL,
    breaks = seq(-2, 14),
    labels = c('<5', '', '10', '', '40', '', '160', '', '640', '', '2560', '', '10240', '', '40960', '', '163840')
  ) +
  coord_cartesian(ylim = c(-3, 14)) +
  scale_color_manual(
    values = animal_colors,
    labels = unlist(lapply(animal_dataset_order, function(i) pretty_labels[i]),
                    use.names=TRUE)
  ) +
  labs(
    x = ''
  )
ggsave('som_figures/fig_s26-28_titer_magnitude/fig_s27_titer_magnitude_by_animal_model.png', plot = gp, width = 22, height = 25)
ggsave('som_figures/fig_s26-28_titer_magnitude/fig_s27_titer_magnitude_by_animal_model.pdf', plot = gp, width = 22, height = 25)


# Make the plot colored by assay
# Re-order the maps, so they are ordered by assay
map_gmts %>%
  mutate(
    map = factor(map, levels = assay_dataset_order),
    assay =  factor(unlist(lapply(as.character(map), function(i) assay[i,]), use.names=TRUE), levels = c('CPE', 'FRNT', 'PV-neut','VSV-PV-neut', 'Microneut', 'PRNT'))
  ) -> map_gmts

plot_adjusted_titer_comparison_gmt(map_gmts, sr_group_gmts, colors = assay_colors, shapes = map_shapes,
                                   included_maps = unique(map_info$map), bar_alpha = 1, ag_order = duplicated_ags,
                                   box_colors = list(
                                     `CPE` = '#ff7f00',
                                     `FRNT` = '#377eb8',
                                     `PV-neut` = '#e41a1c',
                                     `VSV-PV-neut` = '#e3c114',
                                     `PRNT` = '#4daf4a',
                                     `Microneut` = '#984ea3'
                                   ), highlighted_map = c()) -> gp
gp <- gp + theme(
  panel.background = element_rect(fill = 'white', color = 'white'),
  strip.text = element_text(size=25),
  axis.text.x = element_text(size=25),
  axis.text.y = element_text(size=20),
  legend.text=element_text(size=20)
) +
  scale_y_titer(threshold = '<5') +
  scale_y_continuous(
    name = 'Titer',
    minor_breaks = NULL,
    breaks = seq(-2, 14),
    labels = c('<5', '', '10', '', '40', '', '160', '', '640', '', '2560', '', '10240', '', '40960', '', '163840')
  ) +
  coord_cartesian(ylim = c(-3, 14)) +
  scale_color_manual(
    values = assay_colors,
    labels = unlist(lapply(assay_dataset_order, function(i) pretty_labels[i]),
                    use.names=TRUE)
  ) +
  labs(
    x = ''
  )
ggsave('som_figures/fig_s26-28_titer_magnitude/fig_s28_titer_magnitude_by_assay.png', plot = gp, width = 22, height = 25)
ggsave('som_figures/fig_s26-28_titer_magnitude/fig_s28_titer_magnitude_by_assay.pdf', plot = gp, width = 22, height = 25)


