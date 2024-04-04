
rm(list = ls())

library(Racmacs)
library(shape)

map <- read.acmap('data/merged_map/maps/merged_duplicated_antigens_only.ace')

srOutline(map) <- adjustcolor(srOutline(map), alpha.f = 0.5)
map <- rotateMap(map, 4)
map <- translateMap(map, c(1, 0.75))
srOutlineWidth(map) <- 0.5
srSize(map) <- 7

agOutlineWidth(map)[agNames(map) %in% c('D614G', 'B.1.351', 'B.1.617.2', 'BA.1', 'BA.2', 'BA.5')] <- 6

png('som_figures/fig_s53_merged_map_colored_by_variant/fig_s53_merged_map_colored_by_variant.png', width = 7.5, height = 6, units = 'in', res=300,
    pointsize = 18)
par(mar = c(0.3, 0, 0.3, 0))

plot(map, xlim = c(1, 9), ylim = c(-2, 4), fill.alpha = 0.9, plot_labels = FALSE,
     outline.alpha = 0.9, grid.col = '#cfcfcf', plot_stress = FALSE,
     grid.margin.col = 'black', show_error_lines = FALSE, cex=0.75,
     indicate_outliers = FALSE)

coords <- agCoords(map)

lines(c(coords['D614G', 1], coords['B.1.351', 1]), c(coords['D614G', 2], coords['B.1.351', 2]), col = rgb(0, 0, 0, alpha = 0.5), lty = 1, type = 'l', lwd = 2)
lines(c(coords['D614G', 1], coords['B.1.617.2', 1]), c(coords['D614G', 2], coords['B.1.617.2', 2]), col = rgb(0, 0, 0, alpha = 0.5), lty = 1, type = 'l', lwd = 2)
lines(c(coords['D614G', 1], coords['BA.1', 1]), c(coords['D614G', 2], coords['BA.1', 2]), col = rgb(0, 0, 0, alpha = 0.5), lty = 1, type = 'l', lwd = 2)
lines(c(coords['BA.1', 1], coords['BA.2', 1]), c(coords['BA.1', 2], coords['BA.2', 2]), col = rgb(0, 0, 0, alpha = 0.5), lty = 1, type = 'l', lwd = 2)
lines(c(coords['BA.1', 1], coords['BA.5', 1]), c(coords['BA.1', 2], coords['BA.5', 2]), col = rgb(0, 0, 0, alpha = 0.5), lty = 1, type = 'l', lwd = 2)

label_adjustments <- matrix(0, numAntigens(map), 2)
rownames(label_adjustments) <- agNames(map)
label_adjustments['B.1.351',] <- c(-0.6, 0.25)
label_adjustments['P.1',] <- c(-0.7, -0.1)
label_adjustments['B.1.1.7+E484K',] <- c(0.9, -0.05)
label_adjustments['B.1.526+E484K',] <- c(0.95, 0.3)
label_adjustments['B.1.621',] <- c(0.7, 0.1)
label_adjustments['BA.5',] <- c(0.8, 0)
label_adjustments['BA.1',] <- c(-0.2, -0.6)
label_adjustments['BA.2',] <- c(0, 0.6)
label_adjustments['BA.2.12.1',] <- c(0.5, 0.35)
label_adjustments['B.1.617.1',] <- c(0.8, 0)
label_adjustments['BA.1.1',] <- c(0.6, 0)
label_adjustments['B.1.1.7',] <- c(-0.65, 0.1)
label_adjustments['D614G',] <- c(-0.7, 0)
label_adjustments['614D',] <- c(0.6, 0.5)
label_adjustments['B.1.526',] <- c(-0.9, -0.7)
label_adjustments['B.1.526+S477N',] <- c(0.6, 0.6)
label_adjustments['B.1.429',] <- c(-0.7, -0.2)
label_adjustments['B.1.617.2',] <- c(0, -0.65)
label_adjustments['B.1.617.2+K417N',] <- c(1.1, 0)
label_adjustments['C.36.3',] <- c(0.7, 0.1)

labels <- agNames(map)
names(labels) <- agNames(map)
labels['B.1.351'] <- 'B.1.351\nBeta'
labels['P.1'] <- 'P.1\nGamma'
labels['B.1.621'] <- 'B.1.621'
labels['BA.1'] <- 'BA.1\nOmicron'
labels['BA.2'] <- 'BA.2\nOmicron'
labels['B.1.617.1'] <- 'B.1.617.1'
labels['C.37'] <- 'C.37'
labels['B.1.1.7'] <- 'B.1.1.7\nAlpha'
labels['B.1.429'] <- 'B.1.429'
labels['B.1.617.2'] <- 'B.1.617.2\nDelta'
labels['B.1.526+E484K'] <- 'B.1.526+E484K'
labels['BA.5'] <- 'BA.5\nOmicron'
labels['B.1.526+S477N'] <- 'B.1.526+\nS477N'

label_size <- rep(0.7, numAntigens(map))
names(label_size) <- agNames(map)

fonts <- rep(1, numAntigens(map))
names(fonts) <- agNames(map)
bold_font <- 2
fonts['B.1.351'] <- bold_font
fonts['B.1.617.2'] <- bold_font
fonts['D614G'] <- bold_font
fonts['BA.1'] <- bold_font
fonts['BA.2'] <- bold_font
fonts['B.1.1.7'] <- bold_font
fonts['P.1'] <- bold_font
fonts['BA.5'] <- bold_font

text(
  agCoords(map) + label_adjustments,
  cex = label_size,
  label = labels,
  font = fonts
)

dev.off()
