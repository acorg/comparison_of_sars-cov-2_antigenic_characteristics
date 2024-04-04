
# Setup workspace
rm(list = ls())
library(tidyverse)
library(Racmacs)

source('code/plotting/scales.R')
source('code/metadata/common.R')

# Read in data
results <- readRDS('data/merged_map/cross_validation/cross_validation_result.rds')
map <- read.acmap('data/merged_map/maps/merged_duplicated_antigens_only.ace')

# Set detectable results subset
detectable_results <- subset(results, measured_titer_type == 1 & is.finite(residual))

# Plot boxplots split by antigen
agFillScale <- function(map) {
  fills <- agFill(map)
  names(fills) <- agNames(map)
  fills
}

levels(detectable_results$sr_group)[11] <- 'mRNA-1273 vaccinated'
levels(detectable_results$sr_group)[12] <- '3x mRNA-1273 BD01 vaccinated'
levels(detectable_results$sr_group)[13] <- '3x mRNA-1273 BD29 vaccinated'
levels(detectable_results$sr_group)[13] <- '3x mRNA-1273 (6 month) vaccinated'
levels(detectable_results$sr_group)[22] <- 'AZD1222 vaccinated'
levels(detectable_results$sr_group)[24] <- 'AZD1222-BNT162b2 vaccinated'
levels(detectable_results$sr_group)[25] <- 'BNT162b2 vaccinated'

detectable_results %>%
  mutate(
    sr_group = factor(sr_group, levels = c(
      'mRNA-1273 vaccinated', '3x mRNA-1273 BD01 vaccinated', '3x mRNA-1273 BD29 vaccinated',
      '3x mRNA-1273 (6 month) vaccinated', 'AZD1222 vaccinated',
      'AZD1222-BNT162b2 vaccinated', 'BNT162b2 vaccinated', 'D614G convalescent', 'B.1.1.7 convalescent', 'B.1.429 convalescent',
      'B.1.637 convalescent', 'P.1 convalescent', 'B.1.351 convalescent', 'B.1.1.7+E484K convalescent',
      'B.1.621 convalescent', 'B.1.526+E484K convalescent', 'R.1 convalescent', 'P.2 convalescent',
      'B.1+E484K convalescent', 'B.1.617.2 convalescent', 'B.1.617.2+K417N convalescent',
      'B.1.617.1 convalescent', 'C.36.3 convalescent', 'C.37 convalescent', 'BA.1 convalescent',
      'BA.1.1 convalescent', 'BA.2 convalescent')),
      ag_name = factor(ag_name, levels = c('D614G', '614D', 'B.1.1.7', 'B.1.526+S477N', 'B.1.526',
                                           'B.1.429', 'P.1', 'B.1.351', 'B.1.621', 'B.1.1.7+E484K',
                                           'B.1.526+E484K', 'R.1', 'P.2', 'B.1.617.2', 'B.1.617.1',
                                           'B.1.617.2+K417N', 'C.37', 'C.36.3', 'BA.1', 'BA.1.1',
                                           'BA.2.12.1', 'BA.2', 'BA.5'))
  ) %>%
  ggplot(
    aes(
      x = ag_name,
      y = residual,
      fill = ag_name
    )
  ) +
  geom_boxplot(
    lwd = 0.25,
    outlier.shape = NA,
    alpha = 0.8,
    show.legend = FALSE
  ) +
  geom_hline(
    yintercept = 0,
    linetype = 'dashed'
  ) +
  scale_fill_manual(
    values = agFillScale(map)
  ) +
  coord_cartesian(
    y = c(-4, 8)
  ) +
  titerplot_theme() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(
    x = ''
  ) +
  facet_wrap(vars(sr_group))
ggsave('som_figures/fig_s59_measured_vs_predicted_titers_boxplot/fig_s59_measured_vs_predicted_titers_boxplot.png', width = 18, height = 11)
ggsave('som_figures/fig_s59_measured_vs_predicted_titers_boxplot/fig_s59_measured_vs_predicted_titers_boxplot.pdf', width = 18, height = 11)


