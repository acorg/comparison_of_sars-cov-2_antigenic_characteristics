
rm(list = ls())

library(tidyverse)
library(bayestestR)
library(ggpubr)
library(ggrepel)
library(Racmacs)
library(patchwork)

source('code/plotting/map_plotstyles.R')
source('code/plotting/plot_posteriors.R')
source('code/metadata/common.R')


# Read in the animal data
draws_animal <- readRDS('data/titer_analyses/foldchange/slope_calculation_by_animal_assay/slope_calculation_draws_animal.rds')

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
draws_assay <- readRDS('data/titer_analyses/foldchange/slope_calculation_by_animal_assay/slope_calculation_draws_assay.rds')

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
                            rev(c('#d95f02', '#235c3f', '#6383f2')), 'Species slope effect', 
                            xlim = c(-0.01, 1.0)) -> animal_plot 
animal_plot + scale_y_discrete(labels=c('human' = 'Human', 'hamster' = 'Hamster', 'mouse' = 'Mouse')) -> animal_plot


plot_posteriors_hpdi_violin(data_assay, draws_assay_long, rev(c('FRNT', 'LV-PV-neut', 'VSV-PV-neut', 'PRNT', 'Microneut', 'CPE')),
                            rev(c('#003366', '#d95f02', '#CCCC00', '#336600', '#800080', '#008080')), 
                            'Assay slope effect', xlim = c(-0.1, 25)) -> assay_plot

ggarrange(
  animal_plot, assay_plot, ncol = 2, labels = c('A', 'B'), font.label = list(size = 25)
) -> gp_combined
ggsave('som_figures/fig_s41_slope_animal_assay/fig_s41_slope_animal_assay.png', plot = gp_combined,
       width = 10, height = 4, dpi = 300)
ggsave('som_figures/fig_s41_slope_animal_assay/fig_s41_slope_animal_assay.pdf', plot = gp_combined,
       width = 10, height = 4, dpi = 300)

