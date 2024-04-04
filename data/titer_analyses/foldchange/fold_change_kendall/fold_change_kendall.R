
rm(list = ls())

library(tidyverse)
library(Racmacs)
library(DescTools)

source('code/metadata/common.R')
source('code/data_generation/load_maps_for_comparison.R')

## Set the serum groups to look at
sr_groups <- c('mRNA-1273', 'D614G convalescent', 'B.1.1.7 convalescent', 'B.1.351 convalescent',
               'P.1 convalescent',  'B.1.617.2 convalescent', 'BA.1 convalescent')

## Get the fold change table
foldchange_table <- readRDS('data/titer_analyses/foldchange/fold_change_calculation/fold_change_data.rds')

## Get the mean foldchange across all datasets
draws_ag <- readRDS('data/titer_analyses/foldchange/slope_calculation/slope_calculation_ag_folddrops_draws.rds')

maps <- load_maps_for_comparison('', sr_group_standardisation = 'data/metadata/sr_group_standardisation_WT_conv_as_D614G_as614D.csv')
data_orig <- arrange_data(maps, dilution_stepsize = 0)
data <- subsetMap(data_orig, sera = srGroups(data_orig) %in% sr_groups, antigens = agNames(data_orig) %in% duplicated_ags)
srGroups(data) <- factor(srGroups(data), levels = sr_groups)

ag_names_table <- tibble(
  ag_name = rep(agNames(data), numSeraGroups(data)),
  sr_group = rep(levels(srGroups(data)), each = numAntigens(data)),
  variable = sprintf('ag_folddrops[%s,%s]', rep(1:length(agNames(data)), numSeraGroups(data)), rep(1:numSeraGroups(data), each = numAntigens(data)))
)

# Calculate the 95% HPDI
par_ci_ag <- bayestestR::ci(draws_ag, ci = 0.95, method = 'HDI')

# Convert the calculated 95% HPDI to a dataframe for plotting
data_ag <- tibble(variable = par_ci_ag$Parameter, ci_low = par_ci_ag$CI_low, ci_high = par_ci_ag$CI_high, mean_diff = colMeans(draws_ag)[ -c(162:164) ])
data_ag <- left_join(data_ag, ag_names_table)

data_ag %>%
  filter(sr_group %in% sr_groups) -> data_ag


# Calculate correlation coefficients
total <- merge(foldchange_table, data_ag, by=c('sr_group', 'ag_name'))

r <- tibble(map <- character(0), sr_group = character(0), tau = numeric(), ci_lower = numeric(), ci_upper = numeric())
for(map_ in (unique(total$map))) {
  for(sr_g in sr_groups) {
    xxx <- subset(total, map == map_ & sr_group == sr_g)
    if(length(unique(xxx$mean_diff.x)) > 1) {
      t <- KendallTauB(x = xxx$mean_diff.x, y = xxx$mean_diff.y, conf.level = 0.68)
      r <- rbind(r, tibble(map = map_, sr_group = sr_g, tau = t['tau_b'], ci_lower = unname(t['lwr.ci']), ci_upper = unname(t['upr.ci'])))
    }
  }
}

# Save the data
saveRDS(r, 'data/titer_analyses/foldchange/fold_change_kendall/fold_change_kendall.rds')
