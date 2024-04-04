
rm(list = ls())

library(tidyverse)
library(lemon)
library(ggtext)
library(Racmacs)
library(patchwork)

source('code/metadata/common.R')
source('code/plotting/scales.R')
source('code/plotting/plot_fold_change.R')
source('code/data_generation/calculate_fold_change.R')
source('code/data_generation/load_maps_for_comparison.R')
source('code/plotting/map_plotstyles.R')

# Set the serum groups to look at
sr_groups <- c('mRNA-1273', 'D614G convalescent', 'B.1.351 convalescent',
               'B.1.617.2 convalescent', 'BA.1 convalescent')
sr_groups_all <- c('mRNA-1273', 'D614G convalescent', 'B.1.1.7 convalescent', 'B.1.351 convalescent',
                   'P.1 convalescent',  'B.1.617.2 convalescent', 'BA.1 convalescent')

# Get the fold change table
foldchange_table <- readRDS('data/titer_analyses/foldchange/fold_change_calculation/fold_change_data.rds')

# Read in slope data
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
data_eff <- tibble(map = par_ci$Parameter, ci_low = par_ci$CI_low, ci_high = par_ci$CI_high, mean_ = colMeans(draws)[ -c(19:21) ])

# Read in the antigen related data
draws_ag <- readRDS('data/titer_analyses/foldchange/slope_calculation/slope_calculation_ag_folddrops_draws.rds')

# Add information about antigen and serum group
maps <- load_maps_for_comparison('', sr_group_standardisation = 'data/metadata/sr_group_standardisation_WT_conv_as_D614G_as614D.csv')
data_orig <- arrange_data(maps, dilution_stepsize = 0)
data <- subsetMap(data_orig, sera = srGroups(data_orig) %in% sr_groups_all, antigens = agNames(data_orig) %in% duplicated_ags)
srGroups(data) <- factor(srGroups(data), levels = sr_groups_all)

ag_names_table <- tibble(
  ag_name = rep(agNames(data), numSeraGroups(data)),
  sr_group = rep(levels(srGroups(data)), each = numAntigens(data)),
  variable = sprintf('ag_folddrops[%s,%s]', rep(1:length(agNames(data)), numSeraGroups(data)), rep(1:numSeraGroups(data), each = numAntigens(data)))
)

# Calculate the 95% HPDI
par_ci_ag <- bayestestR::ci(draws_ag, ci = 0.95, method = 'HDI')

# Convert the calculated 95% HPDI to a dataframe for plotting
data_ag <- tibble(variable = par_ci_ag$Parameter, ci_low = par_ci_ag$CI_low, ci_high = par_ci_ag$CI_high, mean_ = colMeans(draws_ag)[ -c(162:164) ])
data_ag <- left_join(data_ag, ag_names_table)

data_ag %>%
  filter(sr_group %in% sr_groups) -> data_ag

titrated <- c()
for (m in maps) {
  for (sg in unique(srGroups(m))) {
    titrated <- c(titrated, paste(mapName(m), sg))
  }
}

# Read in the antigen related data
draws_ag_ind <- readRDS('data/titer_analyses/foldchange/slope_calculation/slope_calculation_ag_folddrops_draws.rds')

# Add information about antigen and serum group
ag_names_table_ind <- tibble(
  ag_name = rep(agNames(data), numSeraGroups(data)),
  sr_group = rep(levels(srGroups(data)), each = numAntigens(data)),
  variable = sprintf('ag_folddrops[%s,%s]', rep(1:length(agNames(data)), numSeraGroups(data)), rep(1:numSeraGroups(data), each = numAntigens(data)))
)

# Calculate the 95% HPDI
par_ci_ag_ind <- bayestestR::ci(draws_ag_ind, ci = 0.95, method = 'HDI')

# Convert the calculated 95% HPDI to a dataframe for plotting
data_ag_ind <- tibble(variable = par_ci_ag_ind$Parameter, ci_low = par_ci_ag_ind$CI_low, ci_high = par_ci_ag_ind$CI_high, mean_ = colMeans(draws_ag_ind)[ -c(162:164) ])
data_ag_ind <- left_join(data_ag_ind, ag_names_table_ind)

data_ag_ind %>%
  filter(sr_group %in% sr_groups) -> data_ag_ind


# Make legend
legend_df <- left_join(gather(as_tibble(ag_colors)), gather(as_tibble(ag_shapes)), by='key')
colnames(legend_df) <- c('ag_name', 'color', 'shape')
legend_df$x <- rep(seq(1, 7), each = 5)
legend_df$yy <- rep(seq(5, 1), each = 5)
legend_df %>%
  filter(ag_name %in% duplicated_ags) -> legend_df
legend_df$y <- rev(1:length(duplicated_ags))
names(legend_df$color) <- legend_df$ag_name
names(legend_df$shape) <- legend_df$ag_name

legend_df %>%
  ggplot(
    aes(
      x = x,
      y = yy,
      color = ag_name,
      shape = ag_name
    )
  ) +
  geom_point(
    size = 12
  ) +
  scale_color_manual(
    values = legend_df$color
  ) +
  scale_shape_manual(
    values = legend_df$shape
  ) + theme_void() +
  geom_text(
    aes(
      label = ag_name,
      x = x + 0.1
    ),
    color = 'black',
    hjust = 0,
    size = 11
  ) +
  coord_cartesian(
    xlim = c(-0.5, 7),
    ylim = c(0, 6)
  ) +
  guides(color = 'none', shape = 'none') -> legend_plot


### Plot slope by dataset in black and slope from antigens alone in grey
xloc_foldchange_table <- calculate_xloc_foldchange_table(foldchange_table, data_ag, sr_groups, datasetorder = 'assay')

data_eff_1 <- tibble(map = data_eff$map, ci_low = rep(1, 18), ci_high = rep(1, 18), mean_ = rep(1, 18))

slope_data <- calculate_xloc_sampling(data_ag, data_eff, titrated, sr_groups = sr_groups, datasetorder = 'assay')
slope_data$sr_group = factor(slope_data$sr_group, levels = sr_groups)

slope_data_ag <- calculate_xloc_sampling(data_ag, data_eff_1, titrated, sr_groups = sr_groups, datasetorder = 'assay')
slope_data_ag$sr_group = factor(slope_data_ag$sr_group, levels = sr_groups)

fold_change_by_map_with_slopes_main(xloc_foldchange_table, slope_data, slope_data_ag, c(
  'Emory\n<br><span style = 'font-size:25pt'> </span></br>',
  'Innsbruck\n<br><span style = 'font-size:25pt'> </span></br>',
  'Oxford\n<br><span style = 'font-size:25pt'> </span></br>',
  'Galveston\n<br><span style = 'font-size:25pt'> </span></br>',
  'Madison\n<br><span style = 'font-size:25pt'>(FRNT)</span></br>',
  'WUSTL\n<br><span style = 'font-size:25pt'> </span></br>',
  'Duke\n<br><span style = 'font-size:25pt'> </span></br>',
  'FDA\n<br><span style = 'font-size:25pt'> </span></br>',
  'AMC\n<br><span style = 'font-size:25pt'> </span></br>',
  'EMC<br><span style = 'font-size:25pt'>(VeroE6)</span></br>',
  'EMC<br><span style = 'font-size:25pt'>(Calu-3)</span></br>',
  'Geneva\n<br><span style = 'font-size:25pt'> </span></br>',
  'Charité\n<br><span style = 'font-size:25pt'> </span></br>',
  'EMC<br><span style = 'font-size:25pt'>(PRNT)</span></br>',
  'Mt.<span style = 'font-size:5pt'> </span>Sinai\n<br><span style = 'font-size:25pt'> </span></br>',
  'Madison<br><span style = 'font-size:25pt'>(pooled)</span></br>',
  'Madison<br><span style = 'font-size:25pt'>(unpooled)</span></br>',
  'Maryland<br><span style = 'font-size:25pt'>(pooled)</span></br></span></br>'),
  c(5.99, 8.99, 10.99, 13.99, 14.99)) -> gp

gps <- gp + guides(color = 'none', size = 'none', shape = 'none')
gps + legend_plot + plot_layout(heights = c(2, 0.23), widths = c(0.2, 20), ncol = 1) -> gps_l

ggsave(filename = 'main_figures/fig2_fold_change/fig2_fold_change.png', plot = gps_l, width =  45, height = 30)
ggsave(filename = 'main_figures/fig2_fold_change/fig2_fold_change.pdf', plot = gps_l, width =  45, height = 30)

