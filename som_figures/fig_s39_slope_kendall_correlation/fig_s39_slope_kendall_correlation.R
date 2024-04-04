
rm(list = ls())

library(tidyverse)
library(ggh4x)

source('code/metadata/common.R')

# Read in the data
data <- readRDS('data/titer_analyses/foldchange/fold_change_kendall/fold_change_kendall.rds')

data %>%
  mutate(
    map = unlist(lapply(map, function(i) pretty_labels[i]), use.names=TRUE),
    sr_group = replace(sr_group, sr_group == 'mRNA-1273', 'mRNA-1273/BNT162b2 vaccinated'),
    sr_group = factor(sr_group, levels = rev(c('mRNA-1273/BNT162b2 vaccinated', 'D614G convalescent', 
                                               'B.1.1.7 convalescent', 'B.1.351 convalescent',
                                               'P.1 convalescent',  'B.1.617.2 convalescent', 
                                               'BA.1 convalescent'))),
    negative = case_when(ci_lower < 0  ~ 'negative', 
                         ci_lower > 0.5 ~ 'positive'
    ),
    map = factor(map, unlist(lapply(assay_dataset_order, function(i) pretty_labels[i]), use.names=TRUE))
  ) -> data

data$negative <- data$negative %>% replace_na('intermediate')

design <- '
  ABCDEF
  GHI###
  JK####
  LMN###
  O#####
  PQR###
'

data %>%
  ggplot(
    aes(
      x = tau,
      xmin = ci_lower,
      xmax = ci_upper,
      y = sr_group,
      color = negative
    )
  ) +
  geom_vline(
    xintercept = 0
  ) +
  coord_cartesian(
    xlim = c(-1, 1)
  ) +
  geom_pointrange() +
  theme_bw() +
  theme(
    strip.background = element_rect(colour='white', fill='white'),
    strip.text = element_text(size = 12)
  ) +
  scale_color_manual(
    values = list(
    `positive` = '#fb2a6fff',  # yellow
    `negative` = '#24daffff',  # blue
    `intermediate` = '#9900ffff'  # red
    )
  ) +
  guides(color = 'none') +
  labs(
    y = 'CPE                              Microneut                               PRNT                             VSV-PV-neut                        LV-PV-neut                           FRNT',
    x = "Kendall's tau coefficient"
  ) + facet_manual(vars(map), design = design, axes = 'x') -> gp

gp

ggsave('som_figures/fig_s39_slope_kendall_correlation/fig_s39_slope_kendall_correlation.png', width = 12, height = 10.5)
ggsave('som_figures/fig_s39_slope_kendall_correlation/fig_s39_slope_kendall_correlation.pdf', width = 12, height = 10.5)


