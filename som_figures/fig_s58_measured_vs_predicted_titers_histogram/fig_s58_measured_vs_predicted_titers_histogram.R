
# Setup workspace
rm(list = ls())
library(tidyverse)
library(Racmacs)

source('code/plotting/scales.R')
source('code/metadata/common.R')

results <- readRDS('data/merged_map/cross_validation/cross_validation_result.rds')


# Set detectable results subset
detectable_results <- subset(results, measured_titer_type == 1 & is.finite(residual))

# Plot histogram
detectable_rmse <- sqrt(mean(detectable_results$residual^2))
mean(detectable_results$residual)
sd(detectable_results$residual)

png('som_figures/fig_s58_measured_vs_predicted_titers_histogram/fig_s58_measured_vs_predicted_titers_histogram.png',
    width = 10, height = 8, units = 'in', res=300, pointsize = 18)
hist(detectable_results$residual, 40, freq = F, xlab = 'Measured log titer - Predicted log titer',
       main = sprintf('Detectable titer rmse = %s', round(detectable_rmse, 2)))
dev.off()
