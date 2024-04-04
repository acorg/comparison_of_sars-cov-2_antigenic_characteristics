
rm(list = ls())

library(acutils)
library(tidyverse)
library(ggpubr)
library(bayestestR)
library(titertools)

source('code/metadata/common.R')
source('code/plotting/map_plotstyles.R')
source('code/plotting/plot_titers.R')
source('code/plotting/plot_posteriors.R')
source('code/functions/gmts.R')

# This figure consists of three panels, panel A are the unadjusted titers for the B.1.351 serum group,
# panel B are the modeled dataset magnitude effects, and panel C are the adjusted titers for the B.1.351 serum group.

sr_groups <- c('B.1.351 convalescent')

# Plot all antigens that are in at least two maps, and that have representatives in at least two model animals

## Making panel A

# Read in data
map_info <- readRDS('data/titer_analyses/titer_magnitude/dataset_magnitude_effect/dataset_magnitude_effect_map_info_subset.rds')

# Order maps and change emc map name
map_info %>%
  mutate(
    map = factor(map, levels = animal_dataset_order)
  ) -> map_info

# Calculate the gmts for all serum groups/antigens that should be plotted, per map
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
    animal =  factor(unlist(lapply(as.character(map), function(i) model_animal[i,]), use.names=TRUE), levels = c('human', 'hamster', 'mouse'))
  ) -> map_gmts

# Re-name column headers
names(sr_group_gmts)[names(sr_group_gmts) == 'gmt'] <- 'estimate'
names(sr_group_gmts)[names(sr_group_gmts) == 'gmt_lower'] <- 'estimate_lower'
names(sr_group_gmts)[names(sr_group_gmts) == 'gmt_upper'] <- 'estimate_upper'

# Do the plotting
plot_adjusted_titer_comparison_gmt(map_gmts, sr_group_gmts, colors = animal_colors, shapes = map_shapes,
                                   included_maps = unique(map_info$map), bar_alpha = 1, ag_order = duplicated_ags, 
                                   highlighted_map = c()) -> gp_magnitude

gp_magnitude <- gp_magnitude + theme(
  panel.background = element_rect(fill = 'white', color = 'white'),
  strip.background = element_blank(), strip.text = element_blank(),
  axis.text.x = element_text(size=25),
  axis.text.y = element_text(size=20)
  ) + 
  scale_y_titer(threshold = '<5') +
  coord_cartesian(
    ylim = c(-3, 14)
  ) + 
  scale_y_continuous(
    name = 'Titer',
    minor_breaks = NULL,
    breaks = seq(-2, 14),
    labels = c('<5', '', '10', '', '40', '', '160', '', '640', '', '2560', '', '10240', '', '40960', '', '163840')
  ) +
  labs(x = '') +
  guides(color = 'none')
gp_magnitude
ggsave('main_figures/fig1_titer_magnitude_titer_variability/fig1_titer_magnitude_titer_variability_panelA.png', 
       plot = gp_magnitude,
       width = 20, height = 6.7, dpi = 300)
ggsave('main_figures/fig1_titer_magnitude_titer_variability/fig1_titer_magnitude_titer_variability_panelA.pdf', 
       plot = gp_magnitude,
       width = 20, height = 6.7, dpi = 300)

## Making panel C
# Adjust the titers
map_info %>%
  mutate(
    titer = Racmacs:::reactivity_adjust_titers(
      titers = titer,
      adjustment = -dataset_effects
    )
  ) -> map_info_adjusted

# Calculate the gmts for all serum groups/antigens that should be plotted, per map
map_info_adjusted %>%
  filter(sr_group %in% sr_groups) %>%
  filter(ag_name %in% duplicated_ags) %>%
  group_by(
    map, sr_group, ag_name
  ) %>%
  summarise_gmts() %>%
  filter(!is.na(gmt)) -> map_gmts_adjusted

# Calculate the gmts for all serum groups/antigens that should be plotted across all maps
map_gmts_adjusted %>%
  group_by(
    sr_group, ag_name
  ) %>%
  summarise(gmt=median(gmt), gmt_upper=max(gmt), gmt_lower=min(gmt)) %>%
  ungroup() -> sr_group_gmts_adjusted

# Order the serum groups, add information about animal model
map_gmts_adjusted %>%
  mutate(
    sr_group = factor(sr_group, levels = sr_groups),
    animal =  factor(unlist(lapply(as.character(map), function(i) model_animal[i,]), use.names=TRUE), levels = c('human', 'hamster', 'mouse'))
  ) -> map_gmts_adjusted

# Re-name column headers
names(sr_group_gmts_adjusted)[names(sr_group_gmts_adjusted) == 'gmt'] <- 'estimate'
names(sr_group_gmts_adjusted)[names(sr_group_gmts_adjusted) == 'gmt_lower'] <- 'estimate_lower'
names(sr_group_gmts_adjusted)[names(sr_group_gmts_adjusted) == 'gmt_upper'] <- 'estimate_upper'

# Do the plotting
plot_adjusted_titer_comparison_gmt(map_gmts_adjusted, sr_group_gmts_adjusted, colors = animal_colors, shapes = map_shapes,
                                   included_maps = unique(map_info$map), bar_alpha = 1, ag_order = duplicated_ags,
                                   highlighted_map = c()) -> gp_variability
gp_variability <- gp_variability + theme(
  panel.background = element_rect(fill = 'white', color = 'white'),
  strip.background = element_blank(), strip.text = element_blank(),
  axis.text.x = element_text(size=25),
  axis.text.y = element_text(size=20)) + 
  scale_y_titer(threshold = '<5') +
  scale_y_continuous(
    name = 'Titer',
    minor_breaks = NULL,
    breaks = seq(-2, 14),
    labels = c('<5', '', '10', '', '40', '', '160', '', '640', '', '2560', '', '10240', '', '40960', '', '163840')
  ) +
  coord_cartesian(
    ylim = c(-3, 14)
  ) + labs(x = '') +
  guides(color = 'none')
ggsave('main_figures/fig1_titer_magnitude_titer_variability/fig1_titer_magnitude_titer_variability_panelC.png', 
       plot = gp_variability,
       width = 20, height = 6.7, dpi = 300)
ggsave('main_figures/fig1_titer_magnitude_titer_variability/fig1_titer_magnitude_titer_variability_panelC.pdf', 
       plot = gp_variability,
       width = 20, height = 6.7, dpi = 300)


## Plot panel B

# Read in the data
draws <- readRDS('data/titer_analyses/titer_magnitude/dataset_magnitude_effect/dataset_magnitude_effect_posterior_samples.rds')

# Convert to long format
draws %>%
  pivot_longer(
    cols = c('duke', 'maryland', 'galveston', 'emory',
             'madison_pooled', 'madison_unpooled', 'st_louis',
             'oxford', 'mt_sinai_human', 'emc_prnt',
             'innsbruck', 'charite', 'madison_frnt', 'fda', 'geneva',
             'amc', 'emc_calu', 'emc_vero'),
    names_to = 'map') -> draws_long

# Calculate the 95% HPDI
par_ci <- bayestestR::ci(draws, ci = 0.95, method = 'HDI')

# Convert the calculated 95% HPDI to a dataframe for plotting
data <- tibble(map = par_ci$Parameter, ci_low = par_ci$CI_low, ci_high = par_ci$CI_high, mean_ = colMeans(draws)[ -c(19:21) ])

# Do the plotting
plot_posteriors_hpdi_violin(data, draws_long, rev(animal_dataset_order), assay_colors,
                            xlabel = 'Dataset magnitude effect') -> gp_dataset_effect

gp_dataset_effect_annotated <- gp_dataset_effect + 
  scale_y_discrete(
    limits = rev(animal_dataset_order),
    labels = c('Maryland (pooled)', 'WUSTL',
               'Madison (unpooled)', 'Madison (pooled)', 'Madison (FRNT)',
               'Galveston', 'EMC (Calu-3)', 'EMC (VeroE6)', 'Charité', 'EMC (PRNT)',
               'Geneva', 'AMC', 'Mt. Sinai', 'Oxford', 'Innsbruck', 'FDA',
               'Emory', 'Duke')
  ) +
  geom_hline(
    yintercept = c(2.5, 10.5), 
    color = 'grey40'
  ) + 
  annotate(
    geom = 'text',
    x = -15.5, 
    y = 18, 
    label = 'Human',
    color = 'black',
    size = 5,
    hjust = 0
  ) +
  annotate(
    geom = 'text',
    x = -15.5, 
    y = 10, 
    label = 'Hamster',
    color = 'black',
    size = 5,
    hjust = 0
  ) +
  annotate(
    geom = 'text',
    x = -15.5, 
    y = 2, 
    label = 'Mouse',
    color = 'black',
    size = 5,
    hjust = 0
  ) +
  coord_cartesian(
    xlim = c(-15, 15)
  ) 
ggsave('main_figures/fig1_titer_magnitude_titer_variability/fig1_titer_magnitude_titer_variability_panelB.png', 
       plot = gp_dataset_effect_annotated,
       width = 6, height = 6, dpi = 300)
ggsave('main_figures/fig1_titer_magnitude_titer_variability/fig1_titer_magnitude_titer_variability_panelB.pdf', 
       plot = gp_dataset_effect_annotated,
       width = 6, height = 6, dpi = 300)

# Calculate differences in titer magnitude by animal model
humanMean <- mean(filter(data, map %in% human_maps)$mean_)
hamsterMean <- mean(filter(data, map %in% hamster_maps)$mean_)
mouseMean <- mean(filter(data, map %in% mouse_maps)$mean_)



