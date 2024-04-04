
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
map <- read.acmap('data/individual_datasets/geneva/maps/map.ace')

map <- subsetMap(map, sera = srNames(map)[!srNames(map) %in% c('S13', 'S106')])

srHomologousAgs(map) <- as.list(
  match(homologous_ags_geneva[as.character(srGroups(map))], agNames(map))
)

# Read in serum group gmts
ag_means_individual_effects <- readRDS('data/individual_datasets/geneva/individual_effects/ag_gmt_individual_effects.rds')
sr_group_gmts <- ag_means_individual_effects %>% rename(logtiter = ag_gmt_accounting_for_individual_effect)

# Read in outlier data
outliers <- readRDS('data/individual_datasets/geneva/outliers/outliers.rds')

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
  y_aes = 'logtiter',
  standard_cutoff = FALSE
) -> gp

gpx <- gp + coord_cartesian(
  ylim = c(-5.5, 14.5),
  expand = FALSE
) + scale_y_titer(
  #ymin = -4,
  threshold = '<1',
  logthreshold = -3.321928
)

# Annotate region below detectable
gpx <- gpx +
  annotate(
    'tile',
    x = agNames(map),
    y = -4.321928,
    height = 2,
    fill = 'grey50',
    color = NA,
    alpha = 0.3
  )

# Annotate colors for each antigen
for (n in seq_len(numAntigens(map))) {
  gpx <- gpx +
    annotate(
      'tile',
      x = agNames(map)[n],
      y = -5,
      height = 1,
      fill = agFill(map)[n],
      color = NA
    )
}

ggsave('som_figures/fig_s15_geneva_outlier/fig_s15_geneva_outlier.png', plot = gpx, width = 9, height = 8)
ggsave('som_figures/fig_s15_geneva_outlier/fig_s15_geneva_outlier.pdf', plot = gpx, width = 9, height = 8)






