
rm(list = ls())

library(tidyverse)
library(bayestestR)
library(ggpubr)

source('code/metadata/common.R')
source('code/plotting/map_plotstyles.R')
source('code/plotting/plot_posteriors.R')

## Plot panel A
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
plot_posteriors_hpdi_violin(data, draws_long, rev(assay_dataset_order), animal_colors,
                            xlabel = 'Dataset magnitude effect') -> gp_dataset_effect

gp_dataset_effect_annotated <- gp_dataset_effect + 
  geom_hline(
    yintercept = c(3.5, 4.5, 7.5, 9.5, 12.5), 
    color = 'grey40'
  ) + 
  annotate(
    geom = 'text',
    x = -15.5, 
    y = 18, 
    label = 'FRNT',
    color = '#003366',
    size = 5,
    hjust = 0
  ) +
  annotate(
    geom = 'text',
    x = -15.5, 
    y = 12, 
    label = 'LV-PV-neut',
    color = '#d95f02',
    size = 5,
    hjust = 0
  ) +
  annotate(
    geom = 'text',
    x = -15.5, 
    y = 9, 
    label = 'VSV-PV-neut',
    color = '#CCCC00',
    size = 5,
    hjust = 0
  ) +
  annotate(
    geom = 'text',
    x = -15.5, 
    y = 7, 
    label = 'PRNT',
    color = '#336600',
    size = 5,
    hjust = 0
  )  +
  annotate(
    geom = 'text',
    x = -15.5, 
    y = 4, 
    label = 'Microneut',
    color = '#800080',
    size = 5,
    hjust = 0
  )  +
  annotate(
    geom = 'text',
    x = -15.5, 
    y = 3, 
    label = 'CPE',
    color = '#008080',
    size = 5,
    hjust = 0
  ) +
  coord_cartesian(
    xlim = c(-15, 15)
  ) +
  scale_y_discrete(
    limits = rev(assay_dataset_order),
    labels = c('Maryland', 'Madison (unpooled)', 'Madison (pooled)',
               'Mt. Sinai',
               'EMC (PRNT)', 'Charité', 'Geneva',
               'EMC (Calu-3)', 'EMC (VeroE6)',
               'AMC', 'FDA', 'Duke',
               'Madison (FRNT)', 'WUSTL', 'Galveston',
               'Oxford', 'Innsbruck', 'Emory')
  ) 

ggsave('som_figures/fig_s29_dataset_magnitude_effect/fig_s29_dataset_magnitude_effect_panelA.png', 
       plot = gp_dataset_effect_annotated,
       width = 8, height = 7, dpi = 300)

## Plot panel BC
# Read in the animal data
draws_animal <- readRDS('data/titer_analyses/titer_magnitude/dataset_magnitude_effect_animal_model_assay_differences/dataset_magnitude_effect_animal_model_differences_draws.rds')

# Convert to long format
draws_animal %>%
  pivot_longer(
    cols = c('human', 'hamster', 'mouse'),
    names_to = 'map') %>%
  mutate(
    map = factor(map, levels = rev(c('human', 'hamster', 'mouse')))
  )-> draws_animal_long

# Calculate the 95% HPDI
par_ci_animal <- bayestestR::ci(draws_animal, ci = 0.95, method = 'HDI')

# Convert the calculated 95% HPDI to a dataframe for plotting
data_animal <- tibble(map = par_ci_animal$Parameter, ci_low = par_ci_animal$CI_low, ci_high = par_ci_animal$CI_high, mean_ = colMeans(draws_animal)[ -c(4:6) ])

data_animal %>%
  mutate(
    map = factor(map, levels = rev(c('human', 'hamster', 'mouse')))
  ) -> data_animal

# Read in the assay data
draws_assay <- readRDS('data/titer_analyses/titer_magnitude/dataset_magnitude_effect_animal_model_assay_differences/dataset_magnitude_effect_assay_differences_draws.rds')

colnames(draws_assay) <- c('FRNT', 'LV-PV-neut', 'VSV-PV-neut', 'PRNT', 'Microneut', 'CPE', '.chain', '.iteration', '.draw')

# Convert to long format
draws_assay %>%
  pivot_longer(
    cols = c('FRNT', 'LV-PV-neut', 'VSV-PV-neut', 'PRNT', 'Microneut', 'CPE'),
    names_to = 'map') %>%
  mutate(
    map = factor(map, levels = rev(c('FRNT', 'LV-PV-neut', 'VSV-PV-neut', 'PRNT', 'Microneut', 'CPE')))
  ) -> draws_assay_long

# Calculate the 95% HPDI
par_ci_assay <- bayestestR::ci(draws_assay, ci = 0.95, method = 'HDI')

# Convert the calculated 95% HPDI to a dataframe for plotting
data_assay <- tibble(map = par_ci_assay$Parameter, ci_low = par_ci_assay$CI_low, ci_high = par_ci_assay$CI_high, mean_ = colMeans(draws_assay)[ -c(7:9) ])

data_assay %>%
  mutate(
    map = factor(map, levels = rev(c('FRNT', 'LV-PV-neut', 'VSV-PV-neut', 'PRNT', 'Microneut', 'CPE')))
  ) -> data_assay

plot_posteriors_hpdi_violin(data_animal, draws_animal_long, c('human', 'hamster', 'mouse'), 
                            rev(c('#d95f02', '#235c3f', '#6383f2')), 'Species magnitude effect', 
                            xlim = c(-15, 15)) -> animal_plot
animal_plot + scale_y_discrete(labels=c('human' = 'Human', 'hamster' = 'Hamster', 'mouse' = 'Mouse')) -> animal_plot


plot_posteriors_hpdi_violin(data_assay, draws_assay_long, rev(c('FRNT', 'LV-PV-neut', 'VSV-PV-neut', 'PRNT', 'Microneut', 'CPE')),
                            rev(c('#003366', '#d95f02', '#CCCC00', '#336600', '#800080', '#008080')), 
                            'Assay magnitude effect', xlim = c(-15, 15)) -> assay_plot

ggarrange(
  animal_plot, assay_plot, ncol = 2, labels = c('B', 'C'), font.label = list(size = 25)
) -> gp_combined
ggsave('som_figures/fig_s29_dataset_magnitude_effect/fig_s29_dataset_magnitude_effect_panelBC.png', plot = gp_combined,
       width = 10, height = 4, dpi = 300)


## Plot panel DE
# Read in the animal data
draws_animal <- readRDS('data/titer_analyses/titer_magnitude/dataset_magnitude_effect_animal_model_assay_differences_nt50/estimate_dataset_magnitude_effect_animal_model_assay_differences_nt5_draws_animal.rds')

# Convert to long format
draws_animal %>%
  pivot_longer(
    cols = c('human', 'hamster', 'mouse'),
    names_to = 'map') %>%
  mutate(
    map = factor(map, levels = rev(c('human', 'hamster', 'mouse')))
  )-> draws_animal_long

# Calculate the 95% HPDI
par_ci_animal <- bayestestR::ci(draws_animal, ci = 0.95, method = 'HDI')

# Convert the calculated 95% HPDI to a dataframe for plotting
data_animal <- tibble(map = par_ci_animal$Parameter, ci_low = par_ci_animal$CI_low, ci_high = par_ci_animal$CI_high, mean_ = colMeans(draws_animal)[ -c(4:6) ])

data_animal %>%
  mutate(
    map = factor(map, levels = rev(c('human', 'hamster', 'mouse')))
  ) -> data_animal

# Read in the assay data
draws_assay <- readRDS('data/titer_analyses/titer_magnitude/dataset_magnitude_effect_animal_model_assay_differences_nt50/estimate_dataset_magnitude_effect_animal_model_assay_differences_nt5_draws_assay.rds')

# Change column name
rename(draws_assay, 'LV-PV-neut' = `PV-neut`) -> draws_assay

# Convert to long format
draws_assay %>%
  pivot_longer(
    cols = c('FRNT', 'LV-PV-neut', 'VSV-PV-neut', 'PRNT', 'Microneut', 'CPE'),
    names_to = 'map') %>%
  mutate(
    map = factor(map, levels = rev(c('FRNT', 'LV-PV-neut', 'VSV-PV-neut', 'PRNT', 'Microneut', 'CPE')))
  ) -> draws_assay_long

# Calculate the 95% HPDI
par_ci_assay <- bayestestR::ci(draws_assay, ci = 0.95, method = 'HDI')

# Convert the calculated 95% HPDI to a dataframe for plotting
data_assay <- tibble(map = par_ci_assay$Parameter, ci_low = par_ci_assay$CI_low, ci_high = par_ci_assay$CI_high, mean_ = colMeans(draws_assay)[ -c(7:9) ])

data_assay %>%
  mutate(
    map = factor(map, levels = rev(c('FRNT', 'LV-PV-neut', 'VSV-PV-neut', 'PRNT', 'Microneut', 'CPE')))
  ) -> data_assay

plot_posteriors_hpdi_violin(data_animal, draws_animal_long, c('human', 'hamster', 'mouse'), 
                            rev(c('#d95f02', '#235c3f', '#6383f2')), 'Species magnitude effect', 
                            xlim = c(-15, 15)) -> animal_plot
animal_plot + scale_y_discrete(labels=c('human' = 'Human', 'hamster' = 'Hamster', 'mouse' = 'Mouse')) -> animal_plot


plot_posteriors_hpdi_violin(data_assay, draws_assay_long, rev(c('FRNT', 'LV-PV-neut', 'VSV-PV-neut', 'PRNT', 'Microneut', 'CPE')),
                            rev(c('#003366', '#d95f02', '#CCCC00', '#336600', '#800080', '#008080')), 
                            'Assay magnitude effect', xlim = c(-15, 15)) -> assay_plot

ggarrange(
  animal_plot, assay_plot, ncol = 2, labels = c('D', 'E'), font.label = list(size = 25)
) -> gp_combined
ggsave('som_figures/fig_s29_dataset_magnitude_effect/fig_s29_dataset_magnitude_effect_panelDE.png', plot = gp_combined,
       width = 10, height = 4, dpi = 300)


## Plot panel F
# Read in the titer units data
draws_titer_units <- readRDS('data/titer_analyses/titer_magnitude/dataset_magnitude_effect_titer_units/dataset_magnitude_effect_titer_magnitude_titer_units_titer_units_draws.rds')

# Convert to long format
draws_titer_units %>%
  pivot_longer(
    cols = c('NT50', 'NT90', 'NT99', 'IU/ml'),
    names_to = 'map') %>%
  mutate(
    map = factor(map, levels = rev(c('NT50', 'NT90', 'NT99', 'IU/ml')))
  )-> draws_titer_units_long

# Calculate the 95% HPDI
par_ci_titer_units <- bayestestR::ci(draws_titer_units, ci = 0.95, method = 'HDI')

# Convert the calculated 95% HPDI to a dataframe for plotting
data_titer_units <- tibble(map = par_ci_titer_units$Parameter, ci_low = par_ci_titer_units$CI_low, ci_high = par_ci_titer_units$CI_high, mean_ = colMeans(draws_titer_units)[ -c(5:7) ])

data_titer_units %>%
  mutate(
    map = factor(map, levels = rev(c('NT50', 'NT90', 'NT99', 'IU/ml')))
  ) -> data_titer_units


plot_posteriors_hpdi_violin(data_titer_units, draws_titer_units_long, c('NT50', 'NT90', 'NT99', 'IU/ml'),
                            rev(c('blue1', 'darkorchid1', 'brown1', 'darkgoldenrod1')), 'Titer units effect', 
                            xlim = c(-15, 15)) -> titer_units_plot

ggsave('som_figures/fig_s29_dataset_magnitude_effect/fig_s29_dataset_magnitude_effect_panelF.png', plot = titer_units_plot,
       width = 5, height = 5, dpi = 300)

