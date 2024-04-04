
rm(list = ls())

library(tidyverse)
library(bayestestR)
library(ggpubr)
library(ggrepel)
library(Racmacs)
library(patchwork)

source('code/plotting/map_plotstyles.R')
source('code/plotting/plot_posteriors.R')
source('code/data_generation/load_maps_for_comparison.R')
source('code/functions/map_longinfo.R')
source('code/metadata/common.R')


# Read in the data
draws <- readRDS('data/titer_analyses/foldchange/slope_calculation/slope_calculation_dataset_slope_effect_draws.rds')

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

# Do the plotting by assay
plot_posteriors_hpdi_violin(data, draws_long, rev(assay_dataset_order), animal_colors,
                            xlabel = 'Dataset slope effect',
                            xlim = c(0, 1)) -> gp_dataset_effect_assay

x_text_aln <- 2.01

gp_dataset_effect_assay_annotated <- gp_dataset_effect_assay + 
  scale_y_discrete(
    limits = rev(assay_dataset_order),
    labels = c('Maryland', 'Madison (unpooled)', 'Madison (pooled)',
               'Mt. Sinai',
               'EMC (PRNT)', 'Charité', 'Geneva',
               'EMC (Calu-3)', 'EMC (VeroE6)',
               'AMC', 'FDA', 'Duke',
               'Madison (FRNT)', 'WUSTL', 'Galveston',
               'Oxford', 'Innsbruck', 'Emory')
  ) +
  geom_hline(
    yintercept = c(3.5, 4.5, 7.5, 9.5, 12.5), 
    color = 'grey40'
  ) + 
  annotate(
    geom = 'text',
    x = x_text_aln, 
    y = 18, 
    label = 'FRNT',
    color = '#003366',
    size = 5,
    hjust = 'right'
  ) +
  annotate(
    geom = 'text',
    x = x_text_aln, 
    y = 12, 
    label = 'LV-PV-neut',
    color = '#d95f02',
    size = 5,
    hjust = 'right'
  ) +
  annotate(
    geom = 'text',
    x = x_text_aln, 
    y = 9, 
    label = 'VSV-PV-neut',
    color = '#CCCC00',
    size = 5,
    hjust = 'right'
  ) +
  annotate(
    geom = 'text',
    x = x_text_aln, 
    y = 7, 
    label = 'PRNT',
    color = '#336600',
    size = 5,
    hjust = 'right'
  )  +
  annotate(
    geom = 'text',
    x = x_text_aln, 
    y = 4, 
    label = 'Microneut',
    color = '#800080',
    size = 5,
    hjust = 'right'
  )  +
  annotate(
    geom = 'text',
    x = x_text_aln, 
    y = 3, 
    label = 'CPE',
    color = '#008080',
    size = 5,
    hjust = 'right'
  ) +
  coord_cartesian(
    xlim = c(0, 2)
  )

print(gp_dataset_effect_assay_annotated)

# Do the plotting by animal model
plot_posteriors_hpdi_violin(data, draws_long, rev(animal_dataset_order), assay_colors,
                            xlabel = 'Dataset slope effect',
                            xlim = c(0, 1)) -> gp_dataset_effect_animal

gp_dataset_effect_animal_annotated <- gp_dataset_effect_animal + 
  scale_y_discrete(
    limits = rev(animal_dataset_order),
    labels = c('Maryland', 'WUSTL',
               'Madison (unpooled)', 'Madison (pooled)', 'Madison (FRNT)',
               'Galveston', 'EMC (Calu-3)', 'EMC (VeroE6)', 'Charité', 'EMC (PRNT)',
               'Geneva', 'AMC', 'Mt. Sinai', 'Oxford', 'Innsbruck', 'FDA',
               'Duke', 'Emory')
  ) +
  geom_hline(
    yintercept = c(2.5, 10.5), 
    color = 'grey40'
  ) + 
  annotate(
    geom = 'text',
    x = x_text_aln, 
    y = 18, 
    label = 'Human',
    color = '#d95f02',
    size = 5,
    hjust = 'right'
  ) +
  annotate(
    geom = 'text',
    x = x_text_aln, 
    y = 10, 
    label = 'Hamster',
    color = '#235c3f',
    size = 5,
    hjust = 'right'
  ) +
  annotate(
    geom = 'text',
    x = x_text_aln, 
    y = 2, 
    label = 'Mouse',
    color = '#6383f2',
    size = 5,
    hjust = 'right'
  ) +
  coord_cartesian(
    xlim = c(0, 2)
  )

print(gp_dataset_effect_animal_annotated)

# Plot combined figure for slope effects
ggarrange(
  gp_dataset_effect_assay_annotated, gp_dataset_effect_animal_annotated, ncol = 2, labels = c('A', 'B'), font.label = list(size = 25)
) -> gp_combined
ggsave('som_figures/fig_s40_dataset_slope_effect/fig_s40_dataset_slope_effect.png', plot = gp_combined,
       width = 14, height = 6, dpi = 300)
ggsave('som_figures/fig_s40_dataset_slope_effect/fig_s40_dataset_slope_effect.pdf', plot = gp_combined,
       width = 14, height = 6, dpi = 300)
