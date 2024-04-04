
rm(list = ls())

library(tidyverse)
library(Racmacs)
library(reshape2)
library(ggpubr)

source('code/data_generation/load_maps_for_comparison.R')
source('code/metadata/common.R')

# Read in the maps
maps <- load_maps_for_comparison('', sr_group_standardisation = 'data/metadata/sr_group_standardisation_WT_conv_as_D614G_as614D.csv')

# Function for calculating euclidean distances
euclidean_dist <- function(x, y){
  sqrt(sum((x - y)^2))
}


result <- tibble(
  map = character(0),
  Variant = character(0),
  ratio = numeric(0)
)
for (mapName in list('duke', 'fda', 'innsbruck', 'amc', 'emory', 'geneva',
                     'charite', 'madison_frnt', 'emc_prnt', 'emc_vero', 'emc_calu',
                     'oxford', 'mt_sinai_human', 'st_louis', 'galveston',
                     'maryland', 'madison_unpooled', 'madison_pooled')) {
  map_ <- maps[[mapName]]
  print(mapName(map_))
  # figure out the WT ag (614D vs D614G)
  wt_ag <- unname(wt_ags[mapName(map_)])
  
  # get the coordinates
  coords <- agCoords(map_)
  # calculate the distances
  if ('B.1.351' %in% agNames(map_) & 'B.1.617.2' %in% agNames(map_)) {
    b1351dist <- euclidean_dist(coords[wt_ag, ], coords['B.1.351', ])
    b16172dist <- euclidean_dist(coords[wt_ag, ], coords['B.1.617.2', ])
    
    # calculate the ratios
    b1351ratio <- b16172dist / b1351dist
    
    # store the results
    result <- bind_rows(
      result,
      tibble(
        map = mapName(map_),
        Variant = 'B.1.617.2',
        ratio = b1351ratio
      )
    )
  }
  
  if ('B.1.351' %in% agNames(map_) & 'BA.1' %in% agNames(map_)) {
    ba1dist <- euclidean_dist(coords[wt_ag, ], coords['BA.1', ])
    ba1ratio <- ba1dist / b1351dist
    
    result <- bind_rows(
      result,
      tibble(
        map = mapName(map_),
        Variant = 'BA.1',
        ratio = ba1ratio
      )
    )
  }
  if ('B.1.351' %in% agNames(map_) & 'BA.2' %in% agNames(map_)) {
    ba2dist <- euclidean_dist(coords['BA.1', ], coords['BA.2', ])
    ba2ratio <- ba2dist / b1351dist
    
    result <- bind_rows(
      result,
      tibble(
        map = mapName(map_),
        Variant = 'BA.2',
        ratio = ba2ratio
      )
    )
  }
  
}

result_wide <- dcast(result, map ~ Variant)


# Add assay and animal model information
result$assay <- factor(
  unlist(lapply(as.character(result$map), function(i) assay[i,]),
         use.names=TRUE), levels = c('FRNT', 'LV-PV-neut', 'VSV-PV-neut', 'PRNT', 'Microneut', 'CPE'))
result$model_animal <- factor(
  unlist(lapply(as.character(result$map), function(i) model_animal[i,]),
         use.names=TRUE), levels = c('Human', 'Hamster', 'Mouse'))

result_wide$assay <- factor(
  unlist(lapply(as.character(result_wide$map), function(i) assay[i,]),
         use.names=TRUE), levels = c('FRNT', 'LV-PV-neut', 'VSV-PV-neut', 'PRNT', 'Microneut', 'CPE'))
result_wide$model_animal <- factor(
  unlist(lapply(as.character(result_wide$map), function(i) model_animal[i,]),
         use.names=TRUE), levels = c('Human', 'Hamster', 'Mouse'))

# Plot per dataset
result %>%
  mutate(
    map = factor(map, levels = animal_dataset_order)
  ) %>%
  ggplot(
    aes(
      x = map,
      y = ratio,
      fill = Variant
    )
  ) +
  geom_col(
    position = position_dodge2()
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust=1, vjust = 0.5, size = 12),
    strip.text = element_text(size = 12),
    strip.background = element_rect(colour='white', fill='white')
  ) +
  scale_fill_manual(
    values = list(
      `BA.1` = '#EF3737',
      `BA.2` = '#d10fa2',
      `B.1.617.2` = '#d18652'
    )
  ) +
  guides(fill = 'none') +
  labs(
    y = '(D614G - Variant) /\n(D614G - B.1.351)',
    x = ''
  ) +
  coord_cartesian(
    ylim = c(0, 7)
  ) +
  scale_x_discrete(
    limits = animal_dataset_order,
    labels = unlist(lapply(animal_dataset_order, function(i) pretty_labels[i]),
                    use.names=TRUE)
  ) +
  facet_wrap(~Variant)
ggsave('som_figures/fig_s50_maps_relative_distances/fig_s50_maps_relative_distances_panelA_beta.png', width = 9, height = 3.5)


# Plot per animal model
result %>%
  mutate(
    map = factor(map, levels = animal_dataset_order)
  ) %>%
  ggplot(
    aes(
      x = model_animal,
      y = ratio
    )
  ) +
  geom_boxplot() +
  ggbeeswarm::geom_beeswarm(
    aes(
      color = map,
    ),
    size = 3
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust=1, vjust = 0.5, size = 12),
    strip.text = element_text(size = 12),
    strip.background = element_rect(colour='white', fill='white')
  ) +
  scale_color_manual(
    values = map_colors,
    labels = unlist(lapply(animal_dataset_order, function(i) pretty_labels[i]),
                    use.names=TRUE),
    name = 'Dataset'
  ) +
  coord_cartesian(
    ylim = c(0, 6)
  ) +
  labs(
    y = '(D614G - Variant) /\n(D614G - B.1.351)',
    x = ''
  ) +
  facet_wrap(~Variant) -> animal_model_beta


# Plot per assay
result %>%
  mutate(
    map = factor(map, levels = animal_dataset_order)
  ) %>%
  ggplot(
    aes(
      x = assay,
      y = ratio
    )
  ) +
  geom_boxplot() +
  ggbeeswarm::geom_beeswarm(
    aes(
      color = map,
    ),
    size = 3
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust=1, vjust = 0.5, size = 12),
    strip.text = element_text(size = 12),
    strip.background = element_rect(colour='white', fill='white')
  ) +
  scale_color_manual(
    values = map_colors,
    labels = unlist(lapply(animal_dataset_order, function(i) pretty_labels[i]),
                    use.names=TRUE),
    name = 'Dataset'
  ) +
  coord_cartesian(
    ylim = c(0, 6)
  ) +
  labs(
    y = '(D614G - Variant) /\n(D614G - B.1.351)',
    x = ''
  ) +
  facet_wrap(~Variant) -> assay_beta



# ====== Delta as denominator

result <- tibble(
  map = character(0),
  Variant = character(0),
  ratio = numeric(0)
)
for (mapName in list('duke', 'fda', 'innsbruck', 'amc', 'emory', 'geneva',
                     'charite', 'madison_frnt', 'emc_prnt', 'emc_vero', 'emc_calu',
                     'oxford', 'mt_sinai_human', 'st_louis', 'galveston',
                     'maryland', 'madison_unpooled', 'madison_pooled')) {
  map_ <- maps[[mapName]]
  print(mapName(map_))
  # figure out the WT ag (614D vs D614G)
  wt_ag <- unname(wt_ags[mapName(map_)])
  
  # get the coordinates
  coords <- agCoords(map_)
  # calculate the distances
  if ('B.1.351' %in% agNames(map_) & 'B.1.617.2' %in% agNames(map_)) {
    b1351dist <- euclidean_dist(coords[wt_ag, ], coords['B.1.351', ])
    b16172dist <- euclidean_dist(coords[wt_ag, ], coords['B.1.617.2', ])
    
    # calculate the ratios
    b1351ratio <- b1351dist / b16172dist
    
    # store the results
    result <- bind_rows(
      result,
      tibble(
        map = mapName(map_),
        Variant = 'B.1.351',
        ratio = b1351ratio
      )
    )
  }
  
  if ('B.1.617.2' %in% agNames(map_) & 'BA.1' %in% agNames(map_)) {
    ba1dist <- euclidean_dist(coords[wt_ag, ], coords['BA.1', ])
    ba1ratio <- ba1dist / b16172dist
    
    result <- bind_rows(
      result,
      tibble(
        map = mapName(map_),
        Variant = 'BA.1',
        ratio = ba1ratio
      )
    )
  }
  if ('B.1.617.2' %in% agNames(map_) & 'BA.2' %in% agNames(map_)) {
    ba2dist <- euclidean_dist(coords['BA.1', ], coords['BA.2', ])
    ba2ratio <- ba2dist / b16172dist
    
    result <- bind_rows(
      result,
      tibble(
        map = mapName(map_),
        Variant = 'BA.2',
        ratio = ba2ratio
      )
    )
  }
  
}

result_wide <- dcast(result, map ~ Variant)

# Add assay and animal model information
result$assay <- factor(
  unlist(lapply(as.character(result$map), function(i) assay[i,]),
         use.names=TRUE), levels = c('FRNT', 'LV-PV-neut', 'VSV-PV-neut', 'PRNT', 'Microneut', 'CPE'))
result$model_animal <- factor(
  unlist(lapply(as.character(result$map), function(i) model_animal[i,]),
         use.names=TRUE), levels = c('Human', 'Hamster', 'Mouse'))

result_wide$assay <- factor(
  unlist(lapply(as.character(result_wide$map), function(i) assay[i,]),
         use.names=TRUE), levels = c('FRNT', 'LV-PV-neut', 'VSV-PV-neut', 'PRNT', 'Microneut', 'CPE'))
result_wide$model_animal <- factor(
  unlist(lapply(as.character(result_wide$map), function(i) model_animal[i,]),
         use.names=TRUE), levels = c('Human', 'Hamster', 'Mouse'))

# Plot per dataset
result %>%
  mutate(
    map = factor(map, levels = animal_dataset_order)
  ) %>%
  ggplot(
    aes(
      x = map,
      y = ratio,
      fill = Variant
    )
  ) +
  geom_col(
    position = position_dodge2()
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust=1, vjust = 0.5, size = 12),
    strip.text = element_text(size = 12),
    strip.background = element_rect(colour='white', fill='white')
  ) +
  scale_fill_manual(
    values = list(
      `BA.1` = '#EF3737',
      `BA.2` = '#d10fa2',
      `B.1.351` = '#e7ba52'
    )
  ) +
  guides(fill = 'none') +
  labs(
    y = '(D614G - Variant) /\n(D614G - B.1.617.2)',
    x = ''
  ) +
  scale_x_discrete(
    limits = animal_dataset_order,
    labels = unlist(lapply(animal_dataset_order, function(i) pretty_labels[i]),
                    use.names=TRUE)
  ) +
  coord_cartesian(
    ylim = c(0, 7)
  ) +
  facet_wrap(~Variant)
ggsave('som_figures/fig_s50_maps_relative_distances/fig_s50_maps_relative_distances_panelA_delta.png', width = 9, height = 3.5)

# Plot per animal model
result %>%
  mutate(
    map = factor(map, levels = animal_dataset_order)
  ) %>%
  ggplot(
    aes(
      x = model_animal,
      y = ratio
    )
  ) +
  geom_boxplot() +
  ggbeeswarm::geom_beeswarm(
    aes(
      color = map,
    ),
    size = 3
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust=1, vjust = 0.5, size = 12),
    strip.text = element_text(size = 12),
    strip.background = element_rect(colour='white', fill='white')
  ) +
  scale_color_manual(
    values = map_colors,
    labels = unlist(lapply(animal_dataset_order, function(i) pretty_labels[i]),
                    use.names=TRUE),
    name = 'Dataset'
  ) +
  coord_cartesian(
    ylim = c(0, 7)
  ) +
  labs(
    y = '(D614G - Variant) /\n(D614G - B.1.617.2)',
    x = ''
  ) +
  facet_wrap(~Variant) -> animal_model_delta


# Plot per assay
result %>%
  mutate(
    map = factor(map, levels = animal_dataset_order)
  ) %>%
  ggplot(
    aes(
      x = assay,
      y = ratio
    )
  ) +
  geom_boxplot() +
  ggbeeswarm::geom_beeswarm(
    aes(
      color = map,
    ),
    size = 3
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust=1, vjust = 0.5, size = 12),
    strip.text = element_text(size = 12),
    strip.background = element_rect(colour='white', fill='white')
  ) +
  scale_color_manual(
    values = map_colors,
    labels = unlist(lapply(animal_dataset_order, function(i) pretty_labels[i]),
                    use.names=TRUE),
    name = 'Dataset'
  ) +
  coord_cartesian(
    ylim = c(0, 7)
  ) +
  labs(
    y = '(D614G - Variant) /\n(D614G - B.1.617.2)',
    x = ''
  ) +
  facet_wrap(~Variant) -> assay_delta

# Make composite figures
ggarrange(animal_model_beta, animal_model_delta,
          ncol = 1, nrow = 2, common.legend = TRUE, legend = 'right') -> gp
gp <- gp + theme(panel.background = element_rect(fill='white'))
ggsave('som_figures/fig_s50_maps_relative_distances/fig_s50_maps_relative_distances_panelB.png', plot = gp, width = 10, height = 7)


ggarrange(assay_beta, assay_delta,
          ncol = 1, nrow = 2, common.legend = TRUE, legend = 'right') -> gp
gp <- gp + theme(panel.background = element_rect(fill='white'))
ggsave('som_figures/fig_s50_maps_relative_distances/fig_s50_maps_relative_distances_panelC.png', plot = gp, width = 10, height = 7)

