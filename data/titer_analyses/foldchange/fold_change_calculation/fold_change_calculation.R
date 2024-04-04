
rm(list = ls())

library(Racmacs)
library(tidyverse)
library(titertools)

source('code/metadata/common.R')
source('code/data_generation/calculate_fold_change.R')

# Read in data
map_info <- readRDS('data/titer_analyses/titer_magnitude/dataset_magnitude_effect/dataset_magnitude_effect_map_info.rds')

# Set the serum groups to look at
sr_groups <- c('mRNA-1273', 'D614G convalescent', 'B.1.1.7 convalescent',
               'B.1.351 convalescent', 'P.1 convalescent', 'B.1.617.2 convalescent',
               'BA.1 convalescent')

map_info %>%
  mutate(
    titer = Racmacs:::reactivity_adjust_titers(
      titers = titer,
      adjustment = -dataset_effects
    )
  ) %>%
  subset(sr_group %in% sr_groups) %>%
  subset(ag_name %in% duplicated_ags) %>%
  mutate(
    map = factor(map, levels = animal_dataset_order)
  ) -> map_info_adjusted


homologous_ags <- c(
  'mRNA-1273 duke' = 'D614G',
  'mRNA-1273 maryland' = '614D',
  'mRNA-1273 galveston' = '614D',
  'mRNA-1273 emory' = 'D614G',
  'mRNA-1273 madison_pooled' = 'D614G',
  'mRNA-1273 madison_unpooled' = 'D614G',
  'mRNA-1273 st_louis' = 'D614G',
  'mRNA-1273 oxford' = '614D',
  'mRNA-1273 mt_sinai_mouse' = '614D',
  'mRNA-1273 mt_sinai_human' = '614D',
  'mRNA-1273 emc_prnt' = 'D614G',
  'mRNA-1273 kcl' = '614D',
  'mRNA-1273 innsbruck' = 'D614G',
  'mRNA-1273 fda' = 'D614G',
  'mRNA-1273 geneva' = 'D614G',
  'mRNA-1273 amc' = 'D614G',
  'mRNA-1273 charite' = 'D614G',
  'mRNA-1273 madison_frnt' = '614D',
  'mRNA-1273 emc_calu' = 'D614G',
  'mRNA-1273 emc_vero' = 'D614G',
  'D614G convalescent duke' = 'D614G',
  'D614G convalescent maryland' = '614D',
  'D614G convalescent galveston' = '614D',
  'D614G convalescent emory' = 'D614G',
  'D614G convalescent madison_pooled' = 'D614G',
  'D614G convalescent madison_unpooled' = 'D614G',
  'D614G convalescent st_louis' = 'D614G',
  'D614G convalescent oxford' = '614D',
  'D614G convalescent mt_sinai_mouse' = '614D',
  'D614G convalescent mt_sinai_human' = '614D',
  'D614G convalescent emc_prnt' = 'D614G',
  'D614G convalescent kcl' = '614D',
  'D614G convalescent innsbruck' = 'D614G',
  'D614G convalescent fda' = 'D614G',
  'D614G convalescent geneva' = 'D614G',
  'D614G convalescent amc' = 'D614G',
  'D614G convalescent charite' = 'D614G',
  'D614G convalescent madison_frnt' = '614D',
  'D614G convalescent emc_calu' = 'D614G',
  'D614G convalescent emc_vero' = 'D614G',
  'B.1.351 convalescent duke' = 'B.1.351',
  'B.1.351 convalescent maryland' = 'B.1.351',
  'B.1.351 convalescent galveston' = 'B.1.351',
  'B.1.351 convalescent emory' = 'B.1.351',
  'B.1.351 convalescent madison_pooled' = 'B.1.351',
  'B.1.351 convalescent madison_unpooled' = 'B.1.351',
  'B.1.351 convalescent st_louis' = 'B.1.351',
  'B.1.351 convalescent oxford' = 'B.1.351',
  'B.1.351 convalescent mt_sinai_mouse' = 'B.1.351',
  'B.1.351 convalescent mt_sinai_human' = 'B.1.351',
  'B.1.351 convalescent emc_prnt' = 'B.1.351',
  'B.1.351 convalescent kcl' = 'B.1.351',
  'B.1.351 convalescent innsbruck' = 'B.1.351',
  'B.1.351 convalescent fda' = 'B.1.351',
  'B.1.351 convalescent geneva' = 'B.1.351',
  'B.1.351 convalescent amc' = 'B.1.351',
  'B.1.351 convalescent charite' = 'B.1.351',
  'B.1.351 convalescent madison_frnt' = 'B.1.351',
  'B.1.351 convalescent emc_calu' = 'B.1.351',
  'B.1.351 convalescent emc_vero' = 'B.1.351',
  'B.1.617.2 convalescent duke' = 'B.1.617.2',
  'B.1.617.2 convalescent maryland' = 'B.1.617.2',
  'B.1.617.2 convalescent galveston' = 'B.1.617.2',
  'B.1.617.2 convalescent emory' = 'B.1.617.2',
  'B.1.617.2 convalescent madison_pooled' = 'B.1.617.2',
  'B.1.617.2 convalescent madison_unpooled' = 'B.1.617.2',
  'B.1.617.2 convalescent st_louis' = 'B.1.617.2',
  'B.1.617.2 convalescent oxford' = 'B.1.617.2',
  'B.1.617.2 convalescent mt_sinai_mouse' = 'B.1.617.2',
  'B.1.617.2 convalescent mt_sinai_human' = 'B.1.617.2',
  'B.1.617.2 convalescent emc_prnt' = 'B.1.617.2',
  'B.1.617.2 convalescent kcl' = 'B.1.617.2',
  'B.1.617.2 convalescent innsbruck' = 'B.1.617.2',
  'B.1.617.2 convalescent fda' = 'B.1.617.2',
  'B.1.617.2 convalescent geneva' = 'B.1.617.2',
  'B.1.617.2 convalescent amc' = 'B.1.617.2',
  'B.1.617.2 convalescent charite' = 'B.1.617.2',
  'B.1.617.2 convalescent madison_frnt' = 'B.1.617.2',
  'B.1.617.2 convalescent emc_calu' = 'B.1.617.2',
  'B.1.617.2 convalescent emc_vero' = 'B.1.617.2',
  'P.1 convalescent duke' = 'P.1',
  'P.1 convalescent maryland' = 'P.1',
  'P.1 convalescent galveston' = 'P.1',
  'P.1 convalescent emory' = 'P.1',
  'P.1 convalescent madison_pooled' = 'P.1',
  'P.1 convalescent madison_unpooled' = 'P.1',
  'P.1 convalescent st_louis' = 'P.1',
  'P.1 convalescent oxford' = 'P.1',
  'P.1 convalescent mt_sinai_mouse' = 'P.1',
  'P.1 convalescent mt_sinai_human' = 'P.1',
  'P.1 convalescent emc_prnt' = 'P.1',
  'P.1 convalescent kcl' = 'P.1',
  'P.1 convalescent innsbruck' = 'P.1',
  'P.1 convalescent fda' = 'P.1',
  'P.1 convalescent geneva' = 'P.1',
  'P.1 convalescent amc' = 'P.1',
  'P.1 convalescent charite' = 'P.1',
  'P.1 convalescent madison_frnt' = 'P.1',
  'P.1 convalescent emc_calu' = 'P.1',
  'P.1 convalescent emc_vero' = 'P.1',
  'B.1.1.7 convalescent duke' = 'B.1.1.7',
  'B.1.1.7 convalescent maryland' = 'B.1.1.7',
  'B.1.1.7 convalescent galveston' = 'B.1.1.7',
  'B.1.1.7 convalescent emory' = 'B.1.1.7',
  'B.1.1.7 convalescent madison_pooled' = 'B.1.1.7',
  'B.1.1.7 convalescent madison_unpooled' = 'B.1.1.7',
  'B.1.1.7 convalescent st_louis' = 'B.1.1.7',
  'B.1.1.7 convalescent oxford' = 'B.1.1.7',
  'B.1.1.7 convalescent mt_sinai_mouse' = 'B.1.1.7',
  'B.1.1.7 convalescent mt_sinai_human' = 'B.1.1.7',
  'B.1.1.7 convalescent emc_prnt' = 'B.1.1.7',
  'B.1.1.7 convalescent kcl' = 'B.1.1.7',
  'B.1.1.7 convalescent innsbruck' ='B.1.1.7',
  'B.1.1.7 convalescent fda' = 'B.1.1.7',
  'B.1.1.7 convalescent geneva' = 'B.1.1.7',
  'B.1.1.7 convalescent amc' = 'B.1.1.7',
  'B.1.1.7 convalescent charite' = 'B.1.1.7',
  'B.1.1.7 convalescent madison_frnt' = 'B.1.1.7',
  'B.1.1.7 convalescent emc_calu' = 'B.1.1.7',
  'B.1.1.7 convalescent emc_vero' = 'B.1.1.7',
  'BA.1 convalescent duke' = 'BA.1',
  'BA.1 convalescent maryland' = 'BA.1',
  'BA.1 convalescent galveston' = 'BA.1',
  'BA.1 convalescent emory' = 'BA.1',
  'BA.1 convalescent madison_pooled' = 'BA.1',
  'BA.1 convalescent madison_unpooled' = 'BA.1',
  'BA.1 convalescent st_louis' = 'BA.1',
  'BA.1 convalescent oxford' = 'BA.1',
  'BA.1 convalescent mt_sinai_mouse' = 'BA.1',
  'BA.1 convalescent mt_sinai_human' = 'BA.1',
  'BA.1 convalescent emc_prnt' = 'BA.1',
  'BA.1 convalescent kcl' = 'BA.1',
  'BA.1 convalescent innsbruck' ='BA.1',
  'BA.1 convalescent fda' = 'BA.1',
  'BA.1 convalescent geneva' = 'BA.1',
  'BA.1 convalescent amc' = 'BA.1',
  'BA.1 convalescent charite' = 'BA.1',
  'BA.1 convalescent madison_frnt' = 'BA.1',
  'BA.1 convalescent emc_calu' = 'BA.1',
  'BA.1 convalescent emc_vero' = 'BA.1'
)

# Get the fold change table
foldchange_table <- calculate_fold_change(map_info_adjusted, sr_groups, homologous_ags)

# Save the data
saveRDS(foldchange_table, 'data/titer_analyses/foldchange/fold_change_calculation/fold_change_data.rds')


