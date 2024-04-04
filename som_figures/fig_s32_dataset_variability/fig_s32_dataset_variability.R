
rm(list = ls())

library(tidyverse)
library(bayestestR)
library(ggpubr)

source('code/metadata/common.R')
source('code/plotting/map_plotstyles.R')
source('code/plotting/plot_posteriors.R')

# Read in the data
draws <- readRDS('data/titer_analyses/titer_variability/titer_variability_effect/sr_group_gmts_deviation_sd_draws.rds')

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
x <- point_estimate(draws)

# Convert the calculated 95% HPDI to a dataframe for plotting
data <- tibble(map = par_ci$Parameter, ci_low = par_ci$CI_low, ci_high = par_ci$CI_high, mean_ = colMeans(draws)[ -c(19:21) ], median_ = apply(draws,2,median)[ -c(19:21) ], MAP=x$MAP)

# Do the plotting by assay
plot_posteriors_hpdi_violin(data, draws_long, rev(assay_dataset_order), animal_colors,
                            xlabel = 'Dataset variability effect') -> gp_dataset_effect_assay

gp_dataset_effect_assay_annotated <- gp_dataset_effect_assay + 
  geom_hline(
    yintercept = c(3.5, 4.5, 7.5, 9.5, 12.5), 
    color = 'grey40'
  ) + 
  annotate(
    geom = 'text',
    x = 2.7, 
    y = 18, 
    label = 'FRNT',
    color = '#003366',
    size = 5,
    hjust = 0
  ) +
  annotate(
    geom = 'text',
    x = 2.7, 
    y = 12, 
    label = 'LV-PV-neut',
    color = '#d95f02',
    size = 5,
    hjust = 0
  ) +
  annotate(
    geom = 'text',
    x = 2.7, 
    y = 9, 
    label = 'VSV-PV-neut',
    color = '#CCCC00',
    size = 5,
    hjust = 0
  ) +
  annotate(
    geom = 'text',
    x = 2.7, 
    y = 7, 
    label = 'PRNT',
    color = '#336600',
    size = 5,
    hjust = 0
  )  +
  annotate(
    geom = 'text',
    x = 2.7, 
    y = 4, 
    label = 'Microneut',
    color = '#800080',
    size = 5,
    hjust = 0
  )  +
  annotate(
    geom = 'text',
    x = 2.7, 
    y = 3, 
    label = 'CPE',
    color = '#008080',
    size = 5,
    hjust = 0
  ) +
  coord_cartesian(
    xlim = c(0, 3.5)
  ) + scale_y_discrete(
    limits = rev(assay_dataset_order),
    labels = c('Maryland', 'Madison (unpooled)', 'Madison (pooled)',
               'Mt. Sinai',
               'EMC (PRNT)', 'Charité', 'Geneva',
               'EMC (Calu-3)', 'EMC (VeroE6)',
               'AMC', 'FDA', 'Duke',
               'Madison (FRNT)', 'WUSTL', 'Galveston',
               'Oxford', 'Innsbruck', 'Emory')
  )

print(gp_dataset_effect_assay_annotated)

# Do the plotting by animal model
plot_posteriors_hpdi_violin(data, draws_long, rev(animal_dataset_order), assay_colors,
                            xlabel = 'Dataset variability effect') -> gp_dataset_effect_animal

gp_dataset_effect_animal_annotated <- gp_dataset_effect_animal + 
  geom_hline(
    yintercept = c(2.5, 10.5), 
    color = 'grey40'
  ) + 
  annotate(
    geom = 'text',
    x = 3, 
    y = 18, 
    label = 'Human',
    color = '#d95f02',
    size = 5,
    hjust = 0
  ) +
  annotate(
    geom = 'text',
    x = 3, 
    y = 10, 
    label = 'Hamster',
    color = '#235c3f',
    size = 5,
    hjust = 0
  ) +
  annotate(
    geom = 'text',
    x = 3, 
    y = 2, 
    label = 'Mouse',
    color = '#6383f2',
    size = 5,
    hjust = 0
  ) +
  coord_cartesian(
    xlim = c(0, 3.5)
  ) + scale_y_discrete(
    limits = rev(animal_dataset_order),
    labels = c('Maryland', 'WUSTL',
               'Madison (unpooled)', 'Madison (pooled)',
               'Madison (FRNT)', 'Galveston',
               'EMC (Calu-3)', 'EMC (VeroE6)',
               'Charité', 'EMC (PRNT)',
               'Geneva', 'AMC',  'Mt. Sinai',
               'Oxford', 'Innsbruck', 'FDA',
               'Emory', 'Duke')
  )

print(gp_dataset_effect_animal_annotated)

ggarrange(
  gp_dataset_effect_assay_annotated, gp_dataset_effect_animal_annotated, ncol = 2, labels = c('A', 'B'), font.label = list(size = 25)
) -> gp_combined
ggsave('som_figures/fig_s32_dataset_variability/fig_s32_dataset_variability.png', plot = gp_combined,
       width = 14, height = 6, dpi = 300)
ggsave('som_figures/fig_s32_dataset_variability/fig_s32_dataset_variability.pdf', plot = gp_combined,
       width = 14, height = 6, dpi = 300)


