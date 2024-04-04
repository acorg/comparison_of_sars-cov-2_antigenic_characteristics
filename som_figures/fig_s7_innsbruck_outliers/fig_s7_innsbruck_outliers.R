
# Setup workspace
rm(list = ls())
set.seed(100)
library(tidyverse)
library(Racmacs)
library(patchwork)

# Load functions
source('code/plotting/scales.R')
source('code/functions/map_longinfo.R')
source('code/metadata/sr_group_labels.R')
source('code/metadata/common.R')
source('code/plotting/plot_map_titers.R')
source('code/metadata/homologous_ags.R')

# Read the map
map <- read.acmap('data/individual_datasets/innsbruck/maps/map-OmicronI+II+III-thresholded-single_exposure-P1m1-no-nds.ace')

sr_groups <- c(
  'B.1.617.2 convalescent', 'B.1.617.2 convalescent', 'B.1.617.2 convalescent', 'B.1.617.2 convalescent', 
  'B.1.617.2 convalescent', 'B.1.1.7 convalescent', 'B.1.1.7 convalescent', 'B.1.1.7 convalescent', 
  'B.1.1.7 convalescent', 'B.1.1.7 convalescent', 'B.1.1.7 convalescent', 'B.1.1.7 convalescent', 
  'B.1.1.7 convalescent', 'B.1.1.7 convalescent', 'B.1.351 convalescent', 'B.1.351 convalescent', 
  'B.1.351 convalescent', 'B.1.351 convalescent', 'B.1.351 convalescent', 'B.1.351 convalescent', 
  'mRNA-1273',  'mRNA-1273', 'mRNA-1273', 'mRNA-1273', 'mRNA-1273', 'mRNA-1273', 'mRNA-1273', 'mRNA-1273', 
  'mRNA-1273', 'mRNA-1273', 'AstraZeneca', 'AstraZeneca', 'AstraZeneca', 'AstraZeneca', 'AstraZeneca', 
  'AstraZeneca', 'AstraZeneca', 'AstraZeneca', 'AstraZeneca', 'AstraZeneca', 'AstraZeneca-Pfizer', 
  'AstraZeneca-Pfizer', 'AstraZeneca-Pfizer', 'AstraZeneca-Pfizer', 'AstraZeneca-Pfizer', 
  'AstraZeneca-Pfizer', 'AstraZeneca-Pfizer', 'AstraZeneca-Pfizer', 'AstraZeneca-Pfizer', 
  'AstraZeneca-Pfizer', 'Pfizer', 'Pfizer', 'Pfizer', 'Pfizer', 'Pfizer', 'Pfizer', 'Pfizer', 
  'Pfizer', 'Pfizer', 'Pfizer', 'Pfizer', 'BA.1 convalescent', 'BA.1 convalescent', 'BA.1 convalescent', 
  'BA.1 convalescent', 'BA.1 convalescent', 'BA.1 convalescent', 'BA.1 convalescent', 'BA.1 convalescent', 
  'BA.2 convalescent', 'BA.2 convalescent', 'BA.2 convalescent', 'BA.2 convalescent', 'BA.2 convalescent', 
  'BA.2 convalescent', 'BA.2 convalescent', 'BA.2 convalescent', 'BA.2 convalescent', 'D614G convalescent', 
  'D614G convalescent', 'D614G convalescent', 'D614G convalescent', 'D614G convalescent', 
  'D614G convalescent', 'D614G convalescent', 'D614G convalescent', 'D614G convalescent', 
  'D614G convalescent'
)
srGroups(map) <- factor(sr_groups, levels = unique(sr_groups))

srHomologousAgs(map) <- as.list(
  match(homologous_ags_innsbruck[as.character(srGroups(map))], agNames(map))
)

# Read in serum group gmts
ag_means_individual_effects <- readRDS('data/individual_datasets/innsbruck/individual_effects/ag_gmt_individual_effects.rds')
sr_group_gmts <- ag_means_individual_effects %>% rename(logtiter = ag_gmt_accounting_for_individual_effect)

# Read in outlier data
outliers <- readRDS('data/individual_datasets/innsbruck/outliers/outliers.rds')

# Work out antigen order
ag_order <- sr_group_gmts %>% filter(sr_group == 'D614G convalescent') %>% arrange(-logtiter) %>% select(ag_name) %>% unlist() %>% unname()

plot_map_titers(
  map,
  ag_order = ag_order,
  highlighted_sera = list(
    outliers$higher_than_homologous, 
    outliers$manual
  ),
  highlighted_sera_cols = c('red', 'blue'),
  sr_group_gmts = sr_group_gmts,
  sr_individual_effects = NULL,
  plottype = 'lineplot',
  plot_sr_group_gmts = FALSE,
  y_aes = 'logtiter'
) -> gp

ggsave('som_figures/fig_s7_innsbruck_outliers/fig_s7_innsbruck_outliers.png', plot = gp, width = 9, height = 12)
ggsave('som_figures/fig_s7_innsbruck_outliers/fig_s7_innsbruck_outliers.pdf', plot = gp, width = 9, height = 12)

