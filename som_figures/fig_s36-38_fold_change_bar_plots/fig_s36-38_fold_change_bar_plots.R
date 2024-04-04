
rm(list = ls())

library(Racmacs)
library(tidyverse)
library(ggrepel)
library(reshape2)
library(tidytext)
library(ggpubr)
library(titertools)

source('code/metadata/common.R')
source('code/plotting/plot_fold_change.R')
source('code/functions/gmts.R')
source('code/functions/utils_geoms.R')


# Set the serum groups to look at
sr_groups <- c('mRNA-1273', 'D614G convalescent', 'B.1.1.7 convalescent',
               'B.1.351 convalescent', 'P.1 convalescent',
               'B.1.617.2 convalescent',
               'BA.1 convalescent')

# Get the fold change table
foldchange_table <- readRDS('data/titer_analyses/foldchange/fold_change_calculation/fold_change_data.rds')

foldchange_table %>%
  mutate(
    map = factor(map, levels = c(animal_dataset_order, 'estimated gmt'))
  ) -> foldchange_table

p1 <- make_single_foldchange_plot_no_gmt(subset(foldchange_table, sr_group == 'mRNA-1273'), 
                                         title = 'mRNA-1273 or BNT162b2 vaccinated', cols = map_colors) 
  + scale_color_manual(
    values = map_colors,
    name = '',
    labels = unlist(lapply(animal_dataset_order, function(i) pretty_labels[i]), use.names=TRUE))

p2 <- make_single_foldchange_plot_no_gmt(subset(foldchange_table, sr_group == 'D614G convalescent'), 
                                         title = 'D614G convalescent', cols = map_colors)
  + scale_color_manual(
    values = map_colors,
    name = '',
    labels = unlist(lapply(animal_dataset_order, function(i) pretty_labels[i]), use.names=TRUE))

p3 <- make_single_foldchange_plot_no_gmt(subset(foldchange_table, sr_group == 'B.1.1.7 convalescent'), 
                                         title = 'B.1.1.7 convalescent', cols = map_colors)
  + scale_color_manual(
    values = map_colors,
    name = '',
    labels = unlist(lapply(animal_dataset_order, function(i) pretty_labels[i]), use.names=TRUE))

p4 <- make_single_foldchange_plot_no_gmt(subset(foldchange_table, sr_group == 'B.1.351 convalescent'), 
                                         title = 'B.1.351 convalescent', cols = map_colors)
  + scale_color_manual(
    values = map_colors,
    name = '',
    labels = unlist(lapply(animal_dataset_order, function(i) pretty_labels[i]), use.names=TRUE))

p5 <- make_single_foldchange_plot_no_gmt(subset(foldchange_table, sr_group == 'P.1 convalescent'), 
                                         title = 'P.1 convalescent', cols = map_colors)
  + scale_color_manual(
    values = map_colors,
    name = '',
    labels = unlist(lapply(animal_dataset_order, function(i) pretty_labels[i]), use.names=TRUE))

p6 <- make_single_foldchange_plot_no_gmt(subset(foldchange_table, sr_group == 'B.1.617.2 convalescent'), 
                                         title = 'B.1.617.2 convalescent', cols = map_colors) 
  + scale_color_manual(
    values = map_colors,
    name = '',
    labels = unlist(lapply(animal_dataset_order, function(i) pretty_labels[i]), use.names=TRUE))

p7 <- make_single_foldchange_plot_no_gmt(subset(foldchange_table, sr_group == 'BA.1 convalescent'), 
                                         title = 'BA.1 convalescent', cols = map_colors)
  + scale_color_manual(
    values = map_colors,
    name = '',
    labels = unlist(lapply(animal_dataset_order, function(i) pretty_labels[i]), use.names=TRUE))

ggarrange(p1, p2, p3, p4, p5, p6, p7,
          ncol = 1, nrow = 7, common.legend = TRUE) -> gp
gp <- gp + theme(panel.background = element_rect(fill='white')) 
ggsave('som_figures/fig_s36-38_fold_change_bar_plots/fig_s36_fold_change_bar_plot_by_dataset.png', plot = gp, width = 20, height = 35)
ggsave('som_figures/fig_s36-38_fold_change_bar_plots/fig_s36_fold_change_bar_plot_by_dataset.pdf', plot = gp, width = 20, height = 35)


# Coloured by animal
foldchange_table %>%
  mutate(
    map = factor(map, levels = c(animal_dataset_order, 'estimated gmt'))
  ) -> foldchange_table
p1 <- make_single_foldchange_plot_no_gmt(subset(foldchange_table, sr_group == 'mRNA-1273'), 
                                         title = 'mRNA-1273 or BNT162b2 vaccinated', cols = animal_colors) 
  + scale_color_manual(
    values = animal_colors,
    name = '',
    labels = unlist(lapply(animal_dataset_order, function(i) pretty_labels[i]), use.names=TRUE))

p2 <- make_single_foldchange_plot_no_gmt(subset(foldchange_table, sr_group == 'D614G convalescent'), 
                                         title = 'D614G convalescent', cols = animal_colors) 
  + scale_color_manual(
    values = animal_colors,
    name = '',
    labels = unlist(lapply(animal_dataset_order, function(i) pretty_labels[i]), use.names=TRUE))

p3 <- make_single_foldchange_plot_no_gmt(subset(foldchange_table, sr_group == 'B.1.1.7 convalescent'), 
                                         title = 'B.1.1.7 convalescent', cols = animal_colors) 
  + scale_color_manual(
    values = animal_colors,
    name = '',
    labels = unlist(lapply(animal_dataset_order, function(i) pretty_labels[i]), use.names=TRUE))

p4 <- make_single_foldchange_plot_no_gmt(subset(foldchange_table, sr_group == 'B.1.351 convalescent'), 
                                         title = 'B.1.351 convalescent', cols = animal_colors) 
  + scale_color_manual(
    values = animal_colors,
    name = '',
    labels = unlist(lapply(animal_dataset_order, function(i) pretty_labels[i]), use.names=TRUE))

p5 <- make_single_foldchange_plot_no_gmt(subset(foldchange_table, sr_group == 'P.1 convalescent'), 
                                         title = 'P.1 convalescent', cols = animal_colors) 
  + scale_color_manual(
    values = animal_colors,
    name = '',
    labels = unlist(lapply(animal_dataset_order, function(i) pretty_labels[i]), use.names=TRUE))

p6 <- make_single_foldchange_plot_no_gmt(subset(foldchange_table, sr_group == 'B.1.617.2 convalescent'), 
                                         title = 'B.1.617.2 convalescent', cols = animal_colors) 
  + scale_color_manual(
    values = animal_colors,
    name = '',
    labels = unlist(lapply(animal_dataset_order, function(i) pretty_labels[i]), use.names=TRUE))

p7 <- make_single_foldchange_plot_no_gmt(subset(foldchange_table, sr_group == 'BA.1 convalescent'), 
                                         title = 'BA.1 convalescent', cols = animal_colors) 
  + scale_color_manual(
    values = animal_colors,
    name = '',
    labels = unlist(lapply(animal_dataset_order, function(i) pretty_labels[i]), use.names=TRUE))


ggarrange(p1, p2, p3, p4, p5, p6, p7,
          ncol = 1, nrow = 7, common.legend = TRUE) -> gp
gp <- gp + theme(panel.background = element_rect(fill='white'))
ggsave('som_figures/fig_s36-38_fold_change_bar_plots/fig_s37_fold_change_bar_plot_by_animal_model.png', plot = gp, width = 20, height = 35)
ggsave('som_figures/fig_s36-38_fold_change_bar_plots/fig_s37_fold_change_bar_plot_by_animal_model.pdf', plot = gp, width = 20, height = 35)

# Coloured by assay
foldchange_table %>%
  mutate(
    map = factor(map, levels = c(assay_dataset_order, 'estimated gmt'))
  ) -> foldchange_table
p1 <- make_single_foldchange_plot_no_gmt(subset(foldchange_table, sr_group == 'mRNA-1273'), 
                                         title = 'mRNA-1273 or BNT162b2 vaccinated', 
                                         cols = assay_colors) 
  + scale_color_manual(
    values = assay_colors,
    name = '',
    labels = unlist(lapply(assay_dataset_order, function(i) pretty_labels[i]), use.names=TRUE))

p2 <- make_single_foldchange_plot_no_gmt(subset(foldchange_table, sr_group == 'D614G convalescent'), 
                                         title = 'D614G convalescent', cols = assay_colors) 
  + scale_color_manual(
    values = assay_colors,
    name = '',
    labels = unlist(lapply(assay_dataset_order, function(i) pretty_labels[i]), use.names=TRUE))

p3 <- make_single_foldchange_plot_no_gmt(subset(foldchange_table, sr_group == 'B.1.1.7 convalescent'), 
                                         title = 'B.1.1.7 convalescent', cols = assay_colors) 
  + scale_color_manual(
    values = assay_colors,
    name = '',
    labels = unlist(lapply(assay_dataset_order, function(i) pretty_labels[i]), use.names=TRUE))

p4 <- make_single_foldchange_plot_no_gmt(subset(foldchange_table, sr_group == 'B.1.351 convalescent'), 
                                         title = 'B.1.351 convalescent', cols = assay_colors) 
  + scale_color_manual(
    values = assay_colors,
    name = '',
    labels = unlist(lapply(assay_dataset_order, function(i) pretty_labels[i]), use.names=TRUE))

p5 <- make_single_foldchange_plot_no_gmt(subset(foldchange_table, sr_group == 'P.1 convalescent'), 
                                         title = 'P.1 convalescent', cols = assay_colors) 
  + scale_color_manual(
    values = assay_colors,
    name = '',
    labels = unlist(lapply(assay_dataset_order, function(i) pretty_labels[i]), use.names=TRUE))

p6 <- make_single_foldchange_plot_no_gmt(subset(foldchange_table, sr_group == 'B.1.617.2 convalescent'), 
                                         title = 'B.1.617.2 convalescent', cols = assay_colors) 
  + scale_color_manual(
    values = assay_colors,
    name = '',
    labels = unlist(lapply(assay_dataset_order, function(i) pretty_labels[i]), use.names=TRUE))

p7 <- make_single_foldchange_plot_no_gmt(subset(foldchange_table, sr_group == 'BA.1 convalescent'), 
                                         title = 'BA.1 convalescent', cols = assay_colors) 
  + scale_color_manual(
    values = assay_colors,
    name = '',
    labels = unlist(lapply(assay_dataset_order, function(i) pretty_labels[i]), use.names=TRUE))

ggarrange(p1, p2, p3, p4, p5, p6, p7,
          ncol = 1, nrow = 7, common.legend = T) -> gp
gp <- gp + theme(panel.background = element_rect(fill='white'))
ggsave('som_figures/fig_s36-38_fold_change_bar_plots/fig_s38_fold_change_bar_plot_by_assay.png', plot = gp, width = 20, height = 35)
ggsave('som_figures/fig_s36-38_fold_change_bar_plots/fig_s38_fold_change_bar_plot_by_assay.pdf', plot = gp, width = 20, height = 35)

