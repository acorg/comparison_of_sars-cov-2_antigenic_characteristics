
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
map <- read.acmap('data/individual_datasets/emory/maps/map.ace')

srHomologousAgs(map) <- as.list(
  match(homologous_ags_emory[as.character(srGroups(map))], agNames(map))
)

# Read in serum group gmts
ag_means_individual_effects <- readRDS('data/individual_datasets/emory/individual_effects/ag_gmt_individual_effects.rds')
sr_group_gmts <- ag_means_individual_effects %>% rename(logtiter = ag_gmt_accounting_for_individual_effect)

# Read in outlier data
outliers <- readRDS('data/individual_datasets/emory/outliers/outliers.rds')

# Work out antigen order
ag_order <- sr_group_gmts %>% filter(sr_group == 'WT convalescent') %>% arrange(-logtiter) %>% select(ag_name) %>% unlist() %>% unname()

mapS <- subsetMap(map, sera = srNames(map)[srGroups(map) == 'WT convalescent'])

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

ggsave('som_figures/fig_s3_emory_outliers/fig-s3-emory-outliers.png', plot = gp, width = 9, height = 8)
ggsave('som_figures/fig_s3_emory_outliers/fig-s3-emory-outliers.pdf', plot = gp, width = 9, height = 8)



