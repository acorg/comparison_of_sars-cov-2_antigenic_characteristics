
rm(list = ls())

library(Racmacs)
library(tidyverse)
library(reshape2)
library(ggpubr)
library(titertools)

source('code/metadata/common.R')
source('code/plotting/scales.R')

# Read in data
map_info <- readRDS('data/titer_analyses/titer_magnitude/dataset_magnitude_effect/dataset_magnitude_effect_map_info.rds')

# Adjust titers
map_info %>%
  mutate(
    titer = Racmacs:::reactivity_adjust_titers(
      titers = titer,
      adjustment = -dataset_effects
    )
  ) %>%
  mutate(
    map = factor(map, levels = animal_dataset_order)
  ) -> map_info

sr_groups <- c('mRNA-1273', 'D614G convalescent', 'B.1.1.7 convalescent', 
               'B.1.351 convalescent', 'P.1 convalescent', 
               'B.1.617.2 convalescent', 'BA.1 convalescent')

# Substitutions pairs
subst_pairs <- list(
  c('614D', 'B.1.351'),
  c('D614G', 'B.1.351'),
  c('614D', 'B.1.617.2'),
  c('D614G', 'B.1.617.2'),
  c('614D', 'BA.1'),
  c('D614G', 'BA.1'),
  c('BA.1', 'BA.2'),
  c('BA.2', 'BA.5'),
  c('BA.1', 'BA.5')
)


# Collate data
curated_data <- tibble(
  map = character(0),
  sr_group = character(0),
  ag_a = character(0),
  ag_b = character(0),
  titer_diff = numeric(0),
  titer_diff_upper = numeric(0),
  titer_diff_lower = numeric(0),
  substitution = character(0),
  position = character(0),
  subposition = character(0)
)

for (substitution_pair in subst_pairs) {
  subst1 <- substitution_pair[1]
  subst2 <- substitution_pair[2]
  
  for (m in unique(map_info$map)) {
    for (sg in sr_groups) {

      titers1 <- subset(map_info, map == m & sr_group == sg & ag_name == subst1)$titer
      titers2 <- subset(map_info, map == m & sr_group == sg & ag_name == subst2)$titer
      
      titers <- tibble(titers1, titers2)
      titers <- subset(titers, !(titers1 %in% c('.', '*')) & !(titers2 %in% c('.', '*')))
      
      if (length(titers$titers2) >= 1 & length(titers$titers1) >= 1) {
        result <- titertools::log2diff(
          titers1 = titers$titers1,
          titers2 = titers$titers2,
          dilution_stepsize = 0,
          ci_method = 'HDI',
          mu_prior_mu = 0,
          mu_prior_sigma = 100,
          sigma_prior_alpha = 2,
          sigma_prior_beta = 0.75
        )
        
        if (paste(subst1, subst2) %in% c('614D B.1.351', 'D614G B.1.351', '614D B.1.617.2', 'D614G B.1.617.2')) {
          what = 'wt'
          if (paste(subst1, subst2) %in% c('614D B.1.351', 'D614G B.1.351')) {
            what_sub = 'beta'
          } else {
            what_sub = 'delta'
          }
        }
        if (paste(subst1, subst2) %in% c('614D BA.1', 'D614G BA.1')) {
          what = 'wt-omi'
          what_sub = 'wt-omi'
        }
        if (paste(subst1, subst2) %in% c('BA.1 BA.2', 'BA.2 BA.5', 'BA.1 BA.5')) {
          what = 'omi'
        }
        
        curated_data <- rbind(curated_data, tibble(
          map = m,
          sr_group = sg,
          ag_a = subst1,
          ag_b = subst2,
          titer_diff = result[1],
          titer_diff_upper = result[5],
          titer_diff_lower = result[3],
          substitution = paste(subst1, subst2, sep = ' - '),
          position = what,
          subposition = what_sub
        )
        )                     
      }            
    }
  }
}


# Add column with animal model, assay, and cell type information
curated_data$assay <- factor(
  unlist(lapply(as.character(curated_data$map), function(i) assay[i,]),
         use.names=TRUE), levels = c('FRNT', 'LV-PV-neut', 'VSV-PV-neut', 'PRNT', 'Microneut', 'CPE'))
curated_data$model_animal <- factor(
  unlist(lapply(as.character(curated_data$map), function(i) model_animal[i,]),
         use.names=TRUE), levels = c('Human', 'Hamster', 'Mouse'))
curated_data$cell_type <- factor(
  unlist(lapply(as.character(curated_data$map), function(i) cell_type[i,]),
         use.names=TRUE), levels = c('VeroE6', 'VeroE6-TMPRSS2', 'Vero-TMPRSS2-ACE2', 'HEK293T-TMPRSS2-ACE2', 'HEK293T-ACE2',
                                     'calu-3', 'unknown'))

# Add information about fold drop differences
# Add the distances relative to beta and delta
titer_diff_beta <- c()
titer_diff_delta <- c()
for (i in 1:nrow(curated_data)) {
  row <- curated_data[i, ]
  wt_ag <- unname(wt_ags[as.character(row$map)])
  
  if ('B.1.351' %in% filter(curated_data, map == row$map)$ag_b & 'D614G convalescent' %in% filter(curated_data, map == row$map)$sr_group) {
    beta_diff <- filter(curated_data, map == row$map & sr_group == 'D614G convalescent' & ag_b == 'B.1.351' & ag_a == wt_ag)$titer_diff

    row_diff <- row$titer_diff / beta_diff

    if (length(row_diff) == 0){
      row_diff <- NA
    }

    titer_diff_beta <- c(titer_diff_beta, row_diff)
  } else {
    titer_diff_beta <- c(titer_diff_beta, NA)
  }
  
  
  if ('B.1.617.2' %in% filter(curated_data, map == row$map)$ag_b & 'D614G convalescent' %in% filter(curated_data, map == row$map)$sr_group) {
    delta_diff <- filter(curated_data, map == row$map & sr_group == 'D614G convalescent' & ag_b == 'B.1.617.2' & ag_a == wt_ag)$titer_diff
    
    row_diff <- row$titer_diff / delta_diff
    
    if (length(row_diff) == 0){
      row_diff <- NA
    }
    
    titer_diff_delta <- c(titer_diff_delta, row_diff)
  } else {
    titer_diff_delta <- c(titer_diff_delta, NA)
  }
}

curated_data$titer_diff_beta <- titer_diff_beta
curated_data$titer_diff_delta <- titer_diff_delta

# Order the serum groups
curated_data$sr_group <- factor(curated_data$sr_group, levels = sr_groups)

# Order the maps
curated_data$map = factor(curated_data$map, levels = rev(animal_dataset_order))

# Make box plots
# Subset the data to the relevant comparison (within Omicron, between pre-Omicron and Omicron, and within pre-Omicron)
curated_data %>%
  mutate(
    sr_group = factor(sr_group, levels = sr_groups),
    map = factor(map, levels = animal_dataset_order),
    xx = paste(ag_a, ag_b)
  ) %>%
  filter(!xx %in% c('BA.1 BA.5', 'BA.2 BA.5')) %>%
  filter(sr_group %in% c('BA.1 convalescent')) %>%
  filter(position == 'omi') -> omi

curated_data %>%
  mutate(
    sr_group = factor(sr_group, levels = sr_groups),
    map = factor(map, levels = animal_dataset_order),
    xx = paste(ag_a, ag_b)
  ) %>%
  filter(!xx %in% c('BA.1 BA.5', 'BA.2 BA.5')) %>%
  filter(sr_group %in% c('D614G convalescent')) %>%
  filter(position == 'wt-omi') -> wt_omi

curated_data %>%
  mutate(
    sr_group = factor(sr_group, levels = sr_groups),
    map = factor(map, levels = animal_dataset_order),
    xx = paste(ag_a, ag_b)
  ) %>%
  filter(xx %in% c('D614G B.1.617.2', '614D B.1.617.2')) %>%
  filter(sr_group %in% c('D614G convalescent')) %>%
  filter(subposition == 'delta') -> delta

# Combine the subsetted data
boxplot_subset <- rbind(omi, wt_omi, delta)

# Make the plot by animal model
boxplot_subset %>%
  mutate(
    position = factor(position, levels = c('wt', 'wt-omi', 'omi'))
  ) %>%
  ggplot(
  ) +
  geom_boxplot(
    aes(
      x = model_animal,
      y = titer_diff_beta
    )
  ) +
  geom_point(
    aes(
      x = model_animal,
      y = titer_diff_beta,
      color = map
    ),
    size = 2,
    alpha = 1
  ) + 
  scale_color_manual(
    'Dataset',
    values = map_colors,
    labels = unlist(lapply(animal_dataset_order, function(i) pretty_labels[i]),
                    use.names=TRUE)
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(colour='white', fill='white'),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
  ) +
  labs(
    x = '',
    y = 'Fold change (Variant 1 - Variant 2) /\nFold change (D614G - B.1.351)'
  ) +
  facet_grid(~position, scales = 'free_y', space = 'free_y', 
             labeller = as_labeller(c(`wt` = 'Variant 1: D614G, Variant 2: B.1.617.2,\nSerum group: D614G covalescent', 
                                      `wt-omi` = 'Variant 1: D614G, Variant 2: BA.1,\nSerum group: D614G convalescent', 
                                      `omi` = 'Variant 1: BA.1, Variant 2: BA.2,\nSerum group: BA.1 convalescent'))) -> gp_animal

# Make the plot by assay
boxplot_subset %>%
  mutate(
    position = factor(position, levels = c('wt', 'wt-omi', 'omi'))
  ) %>%
  ggplot(
  ) +
  geom_boxplot(
    aes(
      x = assay,
      y = titer_diff_beta
    )
  ) +
  geom_point(
    aes(
      x = assay,
      y = titer_diff_beta,
      color = map
    ),
    size = 2,
    alpha = 1
  ) + 
  scale_color_manual(
    'Dataset',
    values = map_colors,
    labels = unlist(lapply(animal_dataset_order, function(i) pretty_labels[i]),
                    use.names=TRUE)
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(colour='white', fill='white'),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
  ) +
  labs(
    x = '',
    y = 'Fold change (Variant 1 - Variant 2) /\nFold change (D614G - B.1.351)'
  ) +
  facet_grid(~position, scales = 'free_y', space = 'free_y', 
             labeller = as_labeller(c(`wt` = 'Variant 1: D614G, Variant 2: B.1.617.2,\nSerum group: D614G covalescent', 
                                      `wt-omi` = 'Variant 1: D614G, Variant 2: BA.1,\nSerum group: D614G convalescent', 
                                      `omi` = 'Variant 1: BA.1, Variant 2: BA.2,\nSerum group: BA.1 convalescent'))) -> gp_assay


ggarrange(
  gp_animal, gp_assay, ncol = 1, labels = c('A', 'B'), font.label = list(size = 25), common.legend = TRUE, legend='right'
) -> gp_combined
ggsave('som_figures/fig_s52_relative_fold_change_omicron_pre_omicron/fig_s52_relative_fold_change_omicron_pre_omicron.png', plot = gp_combined,
       width = 8.75, height = 8, dpi = 300)
