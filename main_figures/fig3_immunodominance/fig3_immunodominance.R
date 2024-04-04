
rm(list = ls())

library(Racmacs)
library(tidyverse)
library(reshape2)
library(ggpubr)
library(titertools)
library(ggh4x)

source('code/metadata/common.R')
source('code/plotting/map_plotstyles.R')
source('code/plotting/scales.R')
source('code/functions/gmts.R')

sr_groups <- c('mRNA vaccinated', 'D614G convalescent', 'B.1.1.7 convalescent', 
               'B.1.351 convalescent', 'P.1 convalescent', 
               'B.1.617.2 convalescent')

# Read in data
curated_data <- readRDS('data/titer_analyses/immunodominance/immunodominance_data/immunodominance_data.rds')
levels(curated_data$sr_group)[1] <- 'mRNA vaccinated'

### Plot by dataset with lines by dataset
line_data <- readRDS('data/titer_analyses/immunodominance/immunodominance_data/immunodominance_line_data.rds')
levels(line_data$sr_group)[1] <- 'mRNA vaccinated'

dataset_order <- c('duke', 'fda', 'amc', 'emory', 'innsbruck', 'oxford', 'mt_sinai_human', 'geneva',
                   'galveston', 'madison_frnt', 'emc_prnt', 'charite', 'emc_vero', 'emc_calu',
                   'madison_pooled', 'madison_unpooled', 'st_louis', 'maryland')

line_data %>%
  mutate(
    map = factor(map, levels = dataset_order)
    ) -> line_data_mf

line_data_mf$tmprss2 <- factor(
  unlist(lapply(as.character(line_data_mf$map), function(i) cell_type_2[i,]),
         use.names=TRUE), levels = c('-', 'TMPRSS2'))

## Plot by animal model 484
curated_data %>%
  mutate(
    sr_group = factor(sr_group, levels = c('mRNA vaccinated', 'D614G convalescent', 'B.1.1.7 convalescent', 'B.1.351 convalescent',
                                           'P.1 convalescent',  'B.1.617.2 convalescent', 'BA.1 convalescent')),
    map = factor(map, levels = dataset_order)
  ) %>%
  # Remove data from highly reactive antigens
  filter(!(map == 'madison_pooled' & (substitution %in% c('D614G - B.1.1.7', 'B.1.1.7 - B.1.1.7+E484K', 'B.1.526 - B.1.1.7', '614D - B.1.1.7')))) %>%
  subset(position == '484') %>%
  ggplot(
  ) +
geom_pointrange(
  aes(
    x = titer_diff,
    xmin = titer_diff_lower,
    xmax = titer_diff_upper,
    color = map,
    shape = tmprss2,
    y = substitution
  ),
  size =  0.5,
  alpha = 1,
  position = position_dodge2(
    width = 0.7
  )
) +
  scale_shape_manual(
    values = list(
      `-` = 19,
      `TMPRSS2` = 19
    )
  ) +
  geom_vline(
    xintercept = 0,
    linetype = 'solid',
    color = 'grey40',
    linewidth = 1
  ) +
  geom_vline(
    data = filter(line_data_mf, position == '484' & map %in% unique(filter(curated_data, position == '484')$map)),
    aes(
      xintercept = titer_diff,
      color = map
    ),
    linewidth = 1,
    alpha = 1
  ) +
  coord_cartesian(
    xlim = c(-4, 4)
  ) +
  scale_color_manual(
    'Dataset',
    values = animal_colors,
    labels = unlist(lapply(animal_dataset_order, function(i) pretty_labels[i]),
                    use.names=TRUE)
  ) +
  theme_mapplot() +
  labs(
    x = 'Fold difference (Variant B - Variant A)',
    y = ''
  ) +
  scale_x_continuous(
    breaks = function(x) ceiling(min(x)):floor(max(x)),
    labels = log2foldchange
  ) +
  theme_bw() +
  guides(color='none', shape='none', linetype='none') +
  theme(
    panel.grid.minor.x = element_blank(),
    strip.background = element_rect(colour='white', fill='white'),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.ticks.y=element_blank(),
    axis.text.y=element_blank(),
    strip.text.y=element_blank(),
    strip.text.x=element_text(size = 9),
    legend.text=element_text(size = 11),
    legend.title = element_text(size = 12),
    panel.border = element_rect(color = 'gray70', fill = NA)
  ) +
  facet_grid2(map ~ sr_group, scales = 'free_y', space = 'free_y') -> gp

gp

subset(curated_data, position == '484') -> subset_curated_data
unique(tibble(substitution = subset_curated_data$substitution, map = subset_curated_data$map, `Variant A` = subset_curated_data$ag_a, `Variant B` = subset_curated_data$ag_b)) -> wide_text_data
melt(wide_text_data, id.vars = c('substitution', 'map'))  -> textdata

background_boxes_484 = paste(c('#d95f02', '#d95f02', '#003366', '#003366', '#003366', '#800080', '#336600', '#336600', '#336600', '#008080', '#008080'), '40', sep='')

textdata %>%
  # Remove data from highly reactive antigens
  filter(!(map == 'madison_pooled' & (substitution %in% c('D614G - B.1.1.7', 'B.1.1.7 - B.1.1.7+E484K', 'B.1.526 - B.1.1.7', '614D - B.1.1.7')))) %>%
  mutate(
    map = factor(map, levels = dataset_order)
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
    size = 3.4
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
      size = 10,
      hjust = 0.1,
      face = 'bold'
    ),
    strip.text.y = element_text(
      size = 12,
      face = 'bold',
      vjust = 0.95
    ),
    strip.text.y.left = element_text(angle = 0, hjust=1),
    axis.text.y=element_blank(),
    axis.text.x=element_blank(),
    plot.margin = margin(r = 0, unit = 'cm'),
    panel.spacing.x = unit(-1, 'cm'),
    panel.border = element_blank(),
  ) +
  facet_grid2(map ~ variable, scales = 'free_y', space = 'free_y', switch = 'y', 
              strip = strip_themed(background_y = elem_list_rect(fill = background_boxes_484, colour = background_boxes_484)),
              labeller = labeller(map = pretty_labels)) -> l

ggarrange(l, NULL, gp,
          ncol = 3, nrow = 1, widths = c(0.5, -0.03, 1), align = 'h')

ggsave('main_figures/fig3_immunodominance/fig3_immunodominance_panelA.pdf', width = 13.1, height = 7)


## Plot by animal model 501
curated_data %>%
  mutate(
    sr_group = factor(sr_group, levels = sr_groups),
    map = factor(map, levels = dataset_order)
  ) %>%
  filter(!(map == 'madison_pooled' & (substitution %in% c('D614G - B.1.1.7', 'B.1.1.7 - B.1.1.7+E484K', 'B.1.526 - B.1.1.7', '614D - B.1.1.7')))) %>%
  subset(position == '501') %>%
  ggplot(
  ) +
geom_pointrange(
  aes(
    x = titer_diff,
    xmin = titer_diff_lower,
    xmax = titer_diff_upper,
    color = map,
    shape = tmprss2,
    y = substitution
  ),
  size = 0.5,
  alpha = 1,
  position = position_dodge2(
    width = 0.7
  )
) +
  scale_shape_manual(
    values = list(
      `-` = 19,
      `TMPRSS2` = 19 #21
    )
  ) +
  geom_vline(
    xintercept = 0,
    linetype = 'solid',
    color = 'grey40',
    linewidth = 1
  ) +
  geom_vline(
    data = filter(line_data_mf, position == '501' & map %in% unique(filter(curated_data, position == '501')$map)),
    aes(
      xintercept = titer_diff,
      color = map
    ),
    linewidth = 1,
    alpha = 1
  ) +
  coord_cartesian(
    xlim = c(-4, 4)
  ) +
  scale_color_manual(
    'Dataset',
    values = animal_colors,
    labels = unlist(lapply(animal_dataset_order, function(i) pretty_labels[i]),
                    use.names=TRUE)
  ) +
  theme_mapplot() +
  labs(
    x = 'Fold difference (Variant B - Variant A)',
    y = ''
  ) +
  scale_x_continuous(
    breaks = function(x) ceiling(min(x)):floor(max(x)),
    labels = log2foldchange
  ) +
  theme_bw() +
  guides(color='none', shape='none', linetype='none') +
  theme(
    panel.grid.minor.x = element_blank(),
    strip.background = element_rect(colour='white', fill='white'),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.ticks.y=element_blank(),
    axis.text.y=element_blank(),
    strip.text.y=element_blank(),
    strip.text.x=element_text(size = 9),
    legend.text=element_text(size = 11),
    legend.title = element_text(size = 12),
    panel.border = element_rect(color = 'gray70', fill = NA)
  ) +
  facet_grid2(map ~ sr_group, scales = 'free_y', space = 'free_y') -> gp

gp

subset(curated_data, position == '501') -> subset_curated_data
unique(tibble(substitution = subset_curated_data$substitution, map = subset_curated_data$map, `Variant A` = subset_curated_data$ag_a, `Variant B` = subset_curated_data$ag_b)) -> wide_text_data
melt(wide_text_data, id.vars = c('substitution', 'map'))  -> textdata

background_boxes_501 = paste(c('#d95f02', '#d95f02', '#d95f02', '#003366', '#003366', '#003366', '#800080', '#336600', '#336600', '#336600', '#CCCC00', '#CCCC00', '#008080', '#008080'), '40', sep='')

textdata %>%
  filter(!(map == 'madison_pooled' & (substitution %in% c('D614G - B.1.1.7', 'B.1.1.7 - B.1.1.7+E484K', 'B.1.526 - B.1.1.7', '614D - B.1.1.7')))) %>%
  mutate(
    map = factor(map, levels = dataset_order)
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
    size = 3.4
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
      size = 10,
      hjust = 0.1,
      face = 'bold'
    ),
    strip.text.y = element_text(
      size = 12,
      face = 'bold',
      vjust = 0.95
    ),
    strip.text.y.left = element_text(angle = 0, hjust=1),
    axis.text.y=element_blank(),
    axis.text.x=element_blank(),
    plot.margin = margin(r = 0, unit = 'cm'),
    panel.spacing.x = unit(-1, 'cm'),
    panel.border = element_blank(),
  ) +
  facet_grid2(map ~ variable, scales = 'free_y', space = 'free_y', switch = 'y', 
              strip = strip_themed(background_y = elem_list_rect(fill = background_boxes_501, colour = background_boxes_501)),
              labeller = labeller(map = pretty_labels)) -> l

ggarrange(l, NULL, gp,
          ncol = 3, nrow = 1, widths = c(0.5, -0.03, 1), align = 'h')

ggsave('main_figures/fig3_immunodominance/fig3_immunodominance_panelB.pdf', width = 13.1, height = 7)
