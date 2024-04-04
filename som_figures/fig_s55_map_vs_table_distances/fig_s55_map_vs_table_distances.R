
rm(list = ls())
library(Racmacs)
library(tidyverse)

source('code/functions/map_longinfo.R')
source('code/functions/diagnostics.R')
source('code/functions/imputation.R')
source('code/plotting/scales.R')

# Load the data
map <- read.acmap('data/merged_map/maps/merged_duplicated_antigens_only.ace')

residual_data <- residualErrorTable(map)

# Add antigen and sera info
residual_data %>% 
  left_join(long_ag_info(map), by = 'ag_num') %>%
  left_join(long_sr_info(map), by = 'sr_num') -> residual_data

residual_data$sr_group <- paste(residual_data$sr_group, 'sera', sep=' ')

# Plot the fitted against the measured titers as a scatter plot
residual_data %>%
  ggplot(
    aes(
      x = predicted_logtiter,
      y = measured_logtiter_upper
    )
  ) +
  geom_point(
    alpha = 0.4
  ) +
  theme_bw() +
  scale_x_continuous(
    breaks = -5:14,
    labels = c('<0.625', 2^(-4:14)*10)
  ) +
  scale_y_continuous(
    breaks = -5:14,
    labels = c('<0.625', 2^(-4:14)*10)
  ) +
  coord_cartesian(
    xlim = c(-5, 14),
    ylim = c(-5, 14)
  ) +
  ylab('Measured titers') + 
  xlab('Fitted titers') + 
  geom_abline(intercept = 0, lty=2, col = 'grey') + 
  theme(
    text = element_text(size=15),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) 
ggsave('som_figures/fig_s55_map_vs_table_distances/fig_s55_map_vs_table_distances.png', width = 5.5, height = 5)

