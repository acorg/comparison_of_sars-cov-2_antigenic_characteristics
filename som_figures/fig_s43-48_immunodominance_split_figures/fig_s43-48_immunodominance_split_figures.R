
rm(list = ls())

library(Racmacs)
library(tidyverse)
library(reshape2)
library(ggpubr)
library(titertools)

source('code/metadata/common.R')
source('code/plotting/scales.R')

sr_groups <- c('mRNA-1273 vaccinated', 'D614G convalescent', 'B.1.1.7 convalescent', 
               'B.1.351 convalescent', 'P.1 convalescent', 'B.1.617.2 convalescent')

# Read in data
curated_data <- readRDS('data/titer_analyses/immunodominance/immunodominance_data/immunodominance_data.rds')
levels(curated_data$sr_group)[1] <- 'mRNA vaccinated'

# Add information about assay and animal models
curated_data$assay <- factor(
  unlist(lapply(as.character(curated_data$map), function(i) assay[i,]),
         use.names=TRUE), levels = c('CPE', 'FRNT', 'LV-PV-neut', 'VSV-PV-neut', 'PRNT', 'Microneut'))
curated_data$model_animal <- factor(
  unlist(lapply(as.character(curated_data$map), function(i) model_animal[i,]),
         use.names=TRUE), levels = c('Human', 'Hamster', 'Mouse'))
curated_data$cell_type <- factor(
  unlist(lapply(as.character(curated_data$map), function(i) cell_type[i,]),
         use.names=TRUE), levels = c('VeroE6', 'calu-3', 'HEK293T-ACE2', 'VeroE6-TMPRSS2', 
                                     'Vero-TMPRSS2-ACE2', 'HEK293T-TMPRSS2-ACE2'))

curated_data$map = factor(curated_data$map, levels = rev(assay_dataset_order))

# Plot position 484 split by assay
subset(curated_data, position == '484') %>%
  subset(titer_diff != 0.0) %>%
  ggplot(
    aes(
      y = substitution,
    )
  ) +
  geom_pointrange(
    aes(
      x = titer_diff,
      xmin = titer_diff_lower,
      xmax = titer_diff_upper,
      color = map,
      shape = map
    ),
    size = 0.4,
    position = position_dodge2(
      width = 0.5
    )
  ) +
  coord_cartesian(
    xlim = c(-5, 5)
  ) +
  geom_vline(
    xintercept = 0,
    linetype = 'dashed',
    color = 'grey40'
  ) +
  scale_color_manual(
    'Dataset',
    values = map_colors,
    labels = unlist(lapply(animal_dataset_order, function(i) pretty_labels[i]),
                    use.names=TRUE)
  ) +
  scale_shape_manual(
    'Dataset',
    values = map_shapes_different,
    labels = unlist(lapply(animal_dataset_order, function(i) pretty_labels[i]),
                    use.names=TRUE)
  ) +
  labs(
    x = 'Fold difference (Variant B - Variant A)',
    y = ''
  ) +
  scale_x_continuous(
    breaks = function(x) ceiling(min(x)):floor(max(x)),
    labels = log2foldchange
  ) +
  theme_bw() +
  theme(
    panel.grid.minor.x = element_blank(),
    strip.background = element_rect(colour='white', fill='white'),
    axis.text.y=element_blank(),
    strip.text.y=element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  facet_grid(assay ~ sr_group, margins = 'assay', scales = 'free_x', space = 'free_x') -> gp

subset(curated_data, position == '484') -> subset_curated_data
unique(
  tibble(
    substitution = subset_curated_data$substitution, 
    assay = subset_curated_data$assay, 
    `Variant A` = subset_curated_data$ag_a, 
    `Variant B` = subset_curated_data$ag_b
    )
  ) -> wide_text_data
wide_text_data <- rbind(
  wide_text_data, 
  tibble(
    substitution = '', 
    assay = 'All', 
    `Variant A` = '', 
    `Variant B` = ''
    )
  )
melt(wide_text_data, id.vars = c('substitution', 'assay')) %>%
  complete(nesting(substitution, variable, value), assay) -> textdata

textdata %>%
  filter(
    assay != 'VSV-PV-neut'
    ) %>%
  ggplot(
    aes(
      y = substitution,
      x = 1,
      label = value
    )
  ) +
  geom_text(
    hjust = 0,
    size = 3.1
  ) +
  labs(
    x = '',
    y = ''
  ) +
  coord_cartesian(
    xlim = c(0.98, 1.3),
    ylim = c(2, length(unique(textdata$substitution)))
  ) +
  theme_bw() +
  theme(
    axis.text = element_blank(),
    axis.ticks.y=element_blank(
      
    ),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_text(
      hjust = 0.1,
      face = 'bold'
    ),
    strip.text.y = element_text(
      face = 'bold',
      vjust = 0.95
    ),
    axis.text.y=element_blank(),
    axis.text.x=element_blank(),
    plot.margin = margin(r = 0, unit = 'cm'),
    panel.spacing.x = unit(-1, 'cm'),
    panel.border = element_blank()
  ) +
  facet_grid(assay ~ variable, scales = 'free_y', space = 'free_y', switch = 'y') -> l

ggarrange(l, NULL, gp,
          ncol = 3, nrow = 1, widths = c(0.28, -0.03, 1), align = 'h') -> gpl
gpl
ggsave('som_figures/fig_s43-48_immunodominance_split_figures/fig_s44_immunodominance_484_by_assay.png', plot = gpl, width = 12, height = 12, dpi = 300)
ggsave('som_figures/fig_s43-48_immunodominance_split_figures/fig_s44_immunodominance_484_by_assay.pdf', plot = gpl, width = 12, height = 12, dpi = 300)

# Plot position 484 split by animal model
gp <- gp + facet_grid(model_animal ~ sr_group, margins = 'model_animal', scales = 'free_x', space = 'free_x')

subset(curated_data, position == '484') -> subset_curated_data
unique(
  tibble(
    substitution = subset_curated_data$substitution, 
    model_animal = subset_curated_data$model_animal, 
    `Variant A` = subset_curated_data$ag_a, 
    `Variant B` = subset_curated_data$ag_b
    )
  ) -> wide_text_data
wide_text_data <- rbind(
  wide_text_data, 
  tibble(
    substitution = '', 
    model_animal = 'All', 
    `Variant A` = '', 
    `Variant B` = ''
    )
  )
melt(wide_text_data, id.vars = c('substitution', 'model_animal')) %>%
  complete(nesting(substitution, variable, value), model_animal) -> textdata

textdata %>%
  subset(model_animal != 'Mouse') %>%
  ggplot(
    aes(
      y = substitution,
      x = 1,
      label = value
    )
  ) +
  geom_text(
    hjust = 0,
    size = 3.1
  ) +
  labs(
    x = '',
    y = ''
  ) +
  coord_cartesian(
    xlim = c(0.98, 1.3),
    ylim = c(2, length(unique(textdata$substitution)))
  ) +
  theme_bw() +
  theme(
    axis.text = element_blank(),
    axis.ticks.y=element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_text(
      #size = 10,
      hjust = 0.1,
      face = 'bold'
    ),
    strip.text.y = element_text(
      #size = 15,
      face = 'bold',
      vjust = 0.95
    ),
    #strip.text.y.left = element_text(angle = 0),
    axis.text.y=element_blank(),
    axis.text.x=element_blank(),
    plot.margin = margin(r = 0, unit = 'cm'),
    panel.spacing.x = unit(-1, 'cm'),
    panel.border = element_blank()
  ) +
  facet_grid(model_animal ~ variable, scales = 'free_y', space = 'free_y', switch = 'y') -> l

ggarrange(l, NULL, gp,
          ncol = 3, nrow = 1, widths = c(0.28, -0.03, 1), align = 'h') -> gpl

ggsave('som_figures/fig_s43-48_immunodominance_split_figures/fig_s43_immunodominance_484_by_animal_model.png', plot = gpl, width = 12, height = 6, dpi = 300)
ggsave('som_figures/fig_s43-48_immunodominance_split_figures/fig_s43_immunodominance_484_by_animal_model.pdf', plot = gpl, width = 12, height = 6, dpi = 300)

# Plot position 484 split by cell type
gp <- gp + facet_grid(cell_type ~ sr_group, margins = 'cell_type', scales = 'free_x', space = 'free_x')

subset(curated_data, position == '484') -> subset_curated_data
unique(
  tibble(
    substitution = subset_curated_data$substitution, 
    cell_type = subset_curated_data$cell_type, 
    `Variant A` = subset_curated_data$ag_a, 
    `Variant B` = subset_curated_data$ag_b
    )
  ) -> wide_text_data
wide_text_data <- rbind(wide_text_data, 
                        tibble(
                          substitution = '', 
                          cell_type = 'All', 
                          `Variant A` = '', 
                          `Variant B` = ''
                          )
                        )
melt(wide_text_data, id.vars = c('substitution', 'cell_type')) %>%
  complete(nesting(substitution, variable, value), cell_type) -> textdata

textdata %>%
  subset(cell_type != 'HeLa-ACE2') %>%
  ggplot(
    aes(
      y = substitution,
      x = 1,
      label = value
    )
  ) +
  geom_text(
    hjust = 0,
    size = 3.1
  ) +
  labs(
    x = '',
    y = ''
  ) +
  coord_cartesian(
    xlim = c(0.98, 1.3),
    ylim = c(2, length(unique(textdata$substitution)))
  ) +
  theme_bw() +
  theme(
    axis.text = element_blank(),
    axis.ticks.y=element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_text(
      hjust = 0.1,
      face = 'bold'
    ),
    strip.text.y = element_text(
      face = 'bold',
      vjust = 0.95
    ),
    axis.text.y=element_blank(),
    axis.text.x=element_blank(),
    plot.margin = margin(r = 0, unit = 'cm'),
    panel.spacing.x = unit(-1, 'cm'),
    panel.border = element_blank()
  ) +
  facet_grid(cell_type ~ variable, scales = 'free_y', space = 'free_y', switch = 'y') -> l

ggarrange(l, NULL, gp,
          ncol = 3, nrow = 1, widths = c(0.28, -0.03, 1), align = 'h') -> gpl

ggsave('som_figures/fig_s43-48_immunodominance_split_figures/fig_s47_immunodominance_484_by_cell_type.png', plot = gpl, width = 12, height = 14, dpi = 300)
ggsave('som_figures/fig_s43-48_immunodominance_split_figures/fig_s47_immunodominance_484_by_cell_type.pdf', plot = gpl, width = 12, height = 14, dpi = 300)


# Plot position 501 split by assay
subset(curated_data, position == '501') %>%
  subset(titer_diff != 0.0) %>%
  ggplot(
    aes(
      y = substitution,
    )
  ) +
  geom_pointrange(
    aes(
      x = titer_diff,
      xmin = titer_diff_lower,
      xmax = titer_diff_upper,
      color = map,
      shape = map
    ),
    size = 0.4,
    position = position_dodge2(
      width = 0.5
    )
  ) +
  coord_cartesian(
    xlim = c(-5, 5)
  ) +
  geom_vline(
    xintercept = 0,
    linetype = 'dashed',
    color = 'grey40'
  ) +
  scale_color_manual(
    'Dataset',
    values = map_colors,
    labels = unlist(lapply(animal_dataset_order, function(i) pretty_labels[i]),
                    use.names=TRUE)
  ) +
  scale_shape_manual(
    'Dataset',
    values = map_shapes_different,
    labels = unlist(lapply(animal_dataset_order, function(i) pretty_labels[i]),
                    use.names=TRUE)
  ) +
  labs(
    x = 'Fold difference (Variant B - Variant A)',
    y = ''
  ) +
  scale_x_continuous(
    breaks = function(x) ceiling(min(x)):floor(max(x)),
    labels = log2foldchange
  ) +
  theme_bw() +
  theme(
    panel.grid.minor.x = element_blank(),
    axis.text.y=element_blank(),
    strip.text.y=element_blank(),
    strip.background = element_rect(colour='white', fill='white'),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  facet_grid(assay ~ sr_group, margins = 'assay', scales = 'free_x', space = 'free_x') -> gp


subset(curated_data, position == '501') -> subset_curated_data
unique(
  tibble(
    substitution = subset_curated_data$substitution, 
    assay = subset_curated_data$assay, 
    `Variant A` = subset_curated_data$ag_a, 
    `Variant B` = subset_curated_data$ag_b
    )
  ) -> wide_text_data

wide_text_data <- rbind(wide_text_data, 
                        tibble(substitution = 'D614G - B.1.1.7', 
                               assay = 'All', 
                               `Variant A` = 'D614G', 
                               `Variant B` = 'B.1.1.7'
                               )
                        )
melt(wide_text_data, id.vars = c('substitution', 'assay')) %>%
  complete(nesting(substitution, variable, value), assay) -> textdata

textdata %>%
  ggplot(
    aes(
      y = substitution,
      x = 1,
      label = value
    )
  ) +
  geom_text(
    hjust = 0,
    size = 3.1
  ) +
  labs(
    x = '',
    y = ''
  ) +
  coord_cartesian(
    xlim = c(0.98, 1.3)
  ) +
  theme_bw() +
  theme(
    axis.text = element_blank(),
    axis.ticks.y=element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_text(
      hjust = 0.1,
      face = 'bold'
    ),
    strip.text.y = element_text(
      face = 'bold',
      vjust = 0.95
    ),
    axis.text.y=element_blank(),
    axis.text.x=element_blank(),
    plot.margin = margin(r = 0, unit = 'cm'),
    panel.spacing.x = unit(-0.3, 'cm'),
    panel.border = element_blank()
  ) +
  facet_grid(assay ~ variable, scales = 'free_y', space = 'free_y', switch = 'y') -> l

ggarrange(l, NULL, gp,
          ncol = 3, nrow = 1, widths = c(0.28, -0.02, 1), align = 'h') -> gpl

ggsave('som_figures/fig_s43-48_immunodominance_split_figures/fig_s46_immunodominance_501_by_assay.png', plot = gpl, width = 12, height = 10, dpi = 300)
ggsave('som_figures/fig_s43-48_immunodominance_split_figures/fig_s46_immunodominance_501_by_assay.pdf', plot = gpl, width = 12, height = 10, dpi = 300)

# Plot position 501 split by animal model
gp <- gp + facet_grid(model_animal ~ sr_group, margins = 'model_animal', scales = 'free_x', space = 'free_x')

subset(curated_data, position == '501') -> subset_curated_data
unique(
  tibble(
    substitution = subset_curated_data$substitution, 
    model_animal = subset_curated_data$model_animal, 
    `Variant A` = subset_curated_data$ag_a, 
    `Variant B` = subset_curated_data$ag_b
    )
  ) -> wide_text_data

wide_text_data <- rbind(wide_text_data, 
                        tibble(
                          substitution = 'D614G - B.1.1.7', 
                          model_animal = 'All', 
                          `Variant A` = 'D614G', 
                          `Variant B` = 'B.1.1.7'
                          )
                        )
melt(wide_text_data, id.vars = c('substitution', 'model_animal')) %>%
  complete(nesting(substitution, variable, value), model_animal) -> textdata

textdata %>%
  ggplot(
    aes(
      y = substitution,
      x = 1,
      label = value
    )
  ) +
  geom_text(
    hjust = 0,
    size = 3.1
  ) +
  labs(
    x = '',
    y = ''
  ) +
  coord_cartesian(
    xlim = c(0.98, 1.3)
  ) +
  theme_bw() +
  theme(
    axis.text = element_blank(),
    axis.ticks.y=element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_text(
      hjust = 0.1,
      face = 'bold'
    ),
    strip.text.y = element_text(
      face = 'bold',
      vjust = 0.95
    ),
    axis.text.y=element_blank(),
    axis.text.x=element_blank(),
    plot.margin = margin(r = 0, unit = 'cm'),
    panel.spacing.x = unit(-0.3, 'cm'),
    panel.border = element_blank()
  ) +
  facet_grid(model_animal ~ variable, scales = 'free_y', space = 'free_y', switch = 'y') -> l

ggarrange(l, NULL, gp,
          ncol = 3, nrow = 1, widths = c(0.28, -0.02, 1), align = 'h') -> gpl

ggsave('som_figures/fig_s43-48_immunodominance_split_figures/fig_s45_immunodominance_501_by_animal_model.png', plot = gpl, width = 12, height = 7, dpi = 300)
ggsave('som_figures/fig_s43-48_immunodominance_split_figures/fig_s45_immunodominance_501_by_animal_model.pdf', plot = gpl, width = 12, height = 7, dpi = 300)

# Plot position 501 split by cell type
gp <- gp + facet_grid(cell_type ~ sr_group, margins = 'cell_type', scales = 'free_x', space = 'free_x')

subset(curated_data, position == '501') -> subset_curated_data
unique(
  tibble(
    substitution = subset_curated_data$substitution, 
    cell_type = subset_curated_data$cell_type, 
    `Variant A` = subset_curated_data$ag_a, 
    `Variant B` = subset_curated_data$ag_b
    )
  ) -> wide_text_data

wide_text_data <- rbind(wide_text_data, 
                        tibble(
                          substitution = 'D614G - B.1.1.7', 
                          cell_type = 'All', 
                          `Variant A` = 'D614G', 
                          `Variant B` = 'B.1.1.7'
                          )
                        )
melt(wide_text_data, id.vars = c('substitution', 'cell_type')) %>%
  complete(nesting(substitution, variable, value), cell_type) -> textdata

textdata %>%
  filter(!cell_type %in% c('HeLa-ACE2')) %>%
  ggplot(
    aes(
      y = substitution,
      x = 1,
      label = value
    )
  ) +
  geom_text(
    hjust = 0,
    size = 3.1
  ) +
  labs(
    x = '',
    y = ''
  ) +
  coord_cartesian(
    xlim = c(0.98, 1.3)
  ) +
  theme_bw() +
  theme(
    axis.text = element_blank(),
    axis.ticks.y=element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_text(
      hjust = 0.1,
      face = 'bold'
    ),
    strip.text.y = element_text(
      face = 'bold',
      vjust = 0.95
    ),
    axis.text.y=element_blank(),
    axis.text.x=element_blank(),
    plot.margin = margin(r = 0, unit = 'cm'),
    panel.spacing.x = unit(-0.3, 'cm'),
    panel.border = element_blank()
  ) +
  facet_grid(cell_type ~ variable, scales = 'free_y', space = 'free_y', switch = 'y') -> l

ggarrange(l, NULL, gp,
          ncol = 3, nrow = 1, widths = c(0.28, -0.02, 1), align = 'h') -> gpl

ggsave('som_figures/fig_s43-48_immunodominance_split_figures/fig_s48_immunodominance_501_by_cell_type.png', plot = gpl, width = 12, height = 14, dpi = 300)
ggsave('som_figures/fig_s43-48_immunodominance_split_figures/fig_s48_immunodominance_501_by_cell_type.pdf', plot = gpl, width = 12, height = 14, dpi = 300)

