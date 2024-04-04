
rm(list = ls())

library(Racmacs)
library(tidyverse)
library(reshape2)
library(ggpubr)
library(titertools)

source('code/metadata/common.R')
source('code/plotting/map_plotstyles.R')
source('code/plotting/scales.R')
source('code/functions/gmts.R')
source('code/data_generation/calculate_immunodominance.R')

# Read in data
map_info <- readRDS('data/titer_analyses/titer_magnitude/dataset_magnitude_effect/dataset_magnitude_effect_map_info.rds')

# Adjust titers
map_info %>%
  mutate(
    titer = Racmacs:::reactivity_adjust_titers(
      titers = titer,
      adjustment = -dataset_effects
    )
  ) %>%
  mutate(
    map = factor(map, levels = animal_dataset_order)
  ) -> map_info

sr_groups <- c('mRNA-1273', 'D614G convalescent', 'B.1.1.7 convalescent',
               'B.1.351 convalescent', 'P.1 convalescent',
               'B.1.617.2 convalescent')

# Substitutions pairs
subst_pairs_484 <- list(
  c('614D', 'R.1'),
  c('D614G', 'R.1'),
  c('614D', 'P.2'),
  c('D614G', 'P.2'),
  c('614D', 'B.1.525'),
  c('D614G', 'B.1.525'),
  c('614D', 'B.1.526+E484K'),
  c('D614G', 'B.1.526+E484K'),
  c('D614G', 'B.1+E484K'),
  c('B.1.526', 'R.1'),
  c('B.1.526', 'P.2'),
  c('B.1.526', 'B.1.525'),
  c('B.1.526', 'B.1.526+E484K'),
  c('B.1.1.7', 'B.1.1.7+E484K')
)

subst_pairs_501 <- list(
  c('614D', 'B.1.1.7'),
  c('D614G', 'B.1.1.7'),
  c('B.1.526', 'B.1.1.7'),
  c('B.1.526+E484K', 'B.1.1.7+E484K'),
  c('R.1', 'B.1.1.7+E484K'),
  c('P.2', 'B.1.1.7+E484K'),
  c('B.1.525', 'B.1.1.7+E484K')
)


# Collate data
curated_data <- tibble(
  map = character(0),
  sr_group = character(0),
  ag_a = character(0),
  ag_b = character(0),
  titer_diff = numeric(0),
  titer_diff_upper = numeric(0),
  titer_diff_lower = numeric(0),
  substitution = character(0),
  position = character(0)
)

for (substitution_pair in subst_pairs_484) {
  subst1 <- substitution_pair[1]
  subst2 <- substitution_pair[2]

  for (m in unique(map_info$map)) {
    for (sg in sr_groups) {

      titers1 <- subset(map_info, map == m & sr_group == sg & ag_name == subst1)$titer
      titers2 <- subset(map_info, map == m & sr_group == sg & ag_name == subst2)$titer

      titers <- tibble(titers1, titers2)
      titers <- subset(titers, !(titers1 %in% c('.', '*')) & !(titers2 %in% c('.', '*')))

      if (length(titers$titers2) >= 1 & length(titers$titers1) >= 1) {
        result <- titertools::log2diff(
          titers1 = titers$titers1,
          titers2 = titers$titers2,
          dilution_stepsize = 0,
          ci_method = 'HDI',
          sigma_prior_alpha = 2,
          sigma_prior_beta = 0.75
        )

        curated_data <- rbind(curated_data, tibble(
          map = m,
          sr_group = sg,
          ag_a = subst1,
          ag_b = subst2,
          titer_diff = result[1],
          titer_diff_upper = result[5],
          titer_diff_lower = result[3],
          substitution = paste(subst1, subst2, sep = ' - '),
          position = '484'
        )
        )
      }
    }
  }
}


for (substitution_pair in subst_pairs_501) {
  subst1 <- substitution_pair[1]
  subst2 <- substitution_pair[2]

  for (m in unique(map_info$map)) {
    for (sg in sr_groups) {

      titers1 <- subset(map_info, map == m & sr_group == sg & ag_name == subst1)$titer
      titers2 <- subset(map_info, map == m & sr_group == sg & ag_name == subst2)$titer

      titers <- tibble(titers1, titers2)
      titers <- subset(titers, !(titers1 %in% c('.', '*')) & !(titers2 %in% c('.', '*')))

      if (length(titers$titers2) >= 1 & length(titers$titers1) >= 1) {
        result <- titertools::log2diff(
          titers1 = titers$titers1,
          titers2 = titers$titers2,
          dilution_stepsize = 0,
          ci_method = 'HDI',
          sigma_prior_alpha = 2,
          sigma_prior_beta = 0.75
        )

        curated_data <- rbind(curated_data, tibble(
          map = m,
          sr_group = sg,
          ag_a = subst1,
          ag_b = subst2,
          titer_diff = result[1],
          titer_diff_upper = result[5],
          titer_diff_lower = result[3],
          substitution = paste(subst1, subst2, sep = ' - '),
          position = '501'
        )
        )
      }
    }
  }
}

# Add column with animal model and assay information
curated_data$assay <- factor(
  unlist(lapply(as.character(curated_data$map), function(i) assay[i,]),
         use.names=TRUE), levels = c('CPE', 'FRNT', 'LV-PV-neut', 'VSV-PV-neut', 'PRNT', 'Microneut'))
curated_data$model_animal <- factor(
  unlist(lapply(as.character(curated_data$map), function(i) model_animal[i,]),
         use.names=TRUE), levels = c('Human', 'Hamster', 'Mouse'))
curated_data$tmprss2 <- factor(
  unlist(lapply(as.character(curated_data$map), function(i) cell_type_2[i,]),
         use.names=TRUE), levels = c('-', 'TMPRSS2'))

# Order the serum groups
curated_data$sr_group <- factor(curated_data$sr_group, levels = sr_groups)

# Order the maps
curated_data$map = factor(curated_data$map, levels = animal_dataset_order)

saveRDS(curated_data, 'data/titer_analyses/immunodominance/immunodominance_data/immunodominance_data.rds')


# Summarise fold drop across each map
# Set measurements of B.1.1.7 in madison_pooled to NA to exclude B.1.1.7 due to high reactivity
map_info$titer[map_info$map == 'madison_pooled' & map_info$ag_name == 'B.1.1.7'] <- '*'
map_info$logtiter[map_info$map == 'madison_pooled' & map_info$ag_name == 'B.1.1.7'] <- NA

line_data <- summarise_map_folddrop(map_info, sr_groups)

saveRDS(line_data, 'data/titer_analyses/immunodominance/immunodominance_data/immunodominance_line_data.rds')
