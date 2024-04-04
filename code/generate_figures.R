
## GENERATE MAIN FIGURES

# Generate panels for main figures. Individual panels need to be combined using
# AffinityDesigner. *.afdesign files can be found in the directories containing
# the individual panels for each figure.

message("======== Generating figure 1 (titer magnitude and variability) ========")
source("main_figures/fig1_titer_magnitude_titer_variability/fig1_titer_magnitude_titer_variability.R")

message("======== Generating figure 2 (fold change) ========")
source("main_figures/fig2_fold_change/fig2_fold_change.R")

message("======== Generating figure 3 (immunodominance) ========")
source("main_figures/fig3_immunodominance/fig3_immunodominance.R")

message("======== Generating figure 4 (antigenic map comparison) ========")
source("main_figures/fig4_antigenic_maps/fig4_antigenic_maps.R")


## GENERATE SOM FIGURES

message("======== Generating figures S1, S2, S4, S6, S8, S10, S12, S14, S16-25 (individual antigenic maps) ========")
source("som_figures/fig_s1_2_4_6_8_10_12_14_16-25_individual_map_figures/fig_s1_2_4_6_8_10_12_14_16-25_individual_map_figures.R")

message("======== Generating figure S3 (outliers, emory) ========")
source("som_figures/fig_s3_emory_outliers/fig_s3_emory_outliers.R")

message("======== Generating figure S5 (outliers, fda) ========")
source("som_figures/fig_s5_fda_outliers/fig_s5_fda_outliers.R")

message("======== Generating figure S7 (outliers, innsbruck) ========")
source("som_figures/fig_s7_innsbruck_outliers/fig_s7_innsbruck_outliers.R")

message("======== Generating figure S9 (outliers, oxford) ========")
source("som_figures/fig_s9_oxford_outliers/fig_s9_oxford_outliers.R")

message("======== Generating figure S11 (outliers, mt_sinai) ========")
source("som_figures/fig_s11_mt_sinai_outliers/fig_s11_mt_sinai_outliers.R")

message("======== Generating figure S13 (outliers, amc) ========")
source("som_figures/fig_s13_amc_outliers/fig_s13_amc_outliers.R")

message("======== Generating figure S13 (outliers, geneva) ========")
source("som_figures/fig_s15_geneva_outlier/fig_s15_geneva_outlier.R")

message("======== Generating figures S26-28 (titer magnitude) ========")
source("som_figures/fig_s26-28_titer_magnitude/fig_s26-28_titer_magnitude.R")

message("======== Generating figure S29 (titer magnitude effect) ========")
source("som_figures/fig_s29_dataset_magnitude_effect/fig_s29_dataset_magnitude_effect.R")

message("======== Generating figures S30, S31, S33 (titer variability) ========")
source("som_figures/fig_s30-31-33_titer_variability/fig_s30-31-33_titer_variability.R")

message("======== Generating figure S32 (titer variability) ========")
source("som_figures/fig_s32_dataset_variability/fig_s32_dataset_variability.R")

message("======== Generating figure S34, S35 (fold change) ========")
source("som_figures/fig_s34-35_fold_change/fig_s34-35_fold_change.R")

message("======== Generating figures S36-38 (fold change bar plots) ========")
source("som_figures/fig_s36-38_fold_change_bar_plots/fig_s36-38_fold_change_bar_plots.R")

message("======== Generating figure S39 (fold change, kendall correlation) ========")
source("som_figures/fig_s39_slope_kendall_correlation/fig_s39_slope_kendall_correlation.R")

message("======== Generating figure S40 (dataset slope effect) ========")
source("som_figures/fig_s40_dataset_slope_effect/fig_s40_dataset_slope_effect.R")

message("======== Generating figure S41 (slope effect, animal and assay) ========")
source("som_figures/fig_s41_slope_animal_assay/fig_s41_slope_animal_assay.R")

message("======== Generating figure S42 (immunodominance by assay) ========")
source("som_figures/fig_s42_immunodominance_by_assay/fig_s42_immunodominance_by_assay.R")

message("======== Generating figure S43-48 (immunodominance split figures) ========")
source("som_figures/fig_s43-48_immunodominance_split_figures/fig_s43-48_immunodominance_split_figures.R")

message("======== Generating figure S49 (antigenic maps, panel colored by variant) ========")
source("som_figures/fig_s49_map_panel_colored_by_variant/fig_s49_map_panel_colored_by_variant.R")

message("======== Generating figure S50 (maps, relative distances) ========")
source("som_figures/fig_s50_maps_relative_distances/fig_s50_maps_relative_distances.R")

message("======== Generating figure S51 (antigenic maps, triangulation blobs panel) ========")
source("som_figures/fig_s51_triangulation_blob_panel/fig_s51_triangulation_blob_panel.R")

message("======== Generating figure S52 (relative fold change omicron, pre-omicron) ========")
source("som_figures/fig_s52_relative_fold_change_omicron_pre_omicron/fig_s52_relative_fold_change_omicron_pre_omicron.R")

message("======== Generating figure S53 (merged map colored by variant) ========")
source("som_figures/fig_s53_merged_map_colored_by_variant/fig_s53_merged_map_colored_by_variant.R")

message("======== Generating figure S54 (merged map, dimensionality test) ========")
source("som_figures/fig_s54_merged_map_dimensionality_test/fig_s54_merged_map_dimensionality_test.R")

message("======== Generating figure S55 (merged map, map vs table distances) ========")
source("som_figures/fig_s55_map_vs_table_distances/fig_s55_map_vs_table_distances.R")

message("======== Generating figure S56 (merged map, smooth bootstrap) ========")
source("som_figures/fig_s56_merged_map_smooth_bootstrap/fig_s56_merged_map_smooth_bootstrap.R")

message("======== Generating figure S57 (merged map, resample bootstrap) ========")
source("som_figures/fig_s57_merged_map_resample_bootstrap/fig_s57_merged_map_resample_bootstrap.R")

message("======== Generating figure S58 (merged map, measured vs predicted titers, histogram) ========")
source("som_figures/fig_s58_measured_vs_predicted_titers_histogram/fig_s58_measured_vs_predicted_titers_histogram.R")

message("======== Generating figure S59 (merged map, measured vs predicted titers, boxplot) ========")
source("som_figures/fig_s59_measured_vs_predicted_titers_boxplot/fig_s59_measured_vs_predicted_titers_boxplot.R")

message("======== Generating figure S60 (merged map vs individual maps cross validation) ========")
source("som_figures/fig_s60_merged_map_individual_map_cross_validation/fig_s60_merged_map_individual_map_cross_validation.R")

message("======== Generating figure S61 (merged map vs cpe cross validation) ========")
source("som_figures/fig_s61_merged_map_cpe_maps_cross_validation/fig_s61_merged_map_cpe_maps_cross_validation.R")

message("======== Generating figure S62 (merged map vs individual maps, procrustes) ========")
source("som_figures/fig_s62_merged_map_individual_maps_procrustes/fig_s62_merged_map_individual_maps_procrustes.R")

message("======== Generating figure S63 (individual maps vs merged map, procrustes) ========")
source("som_figures/fig_s63_individual_maps_merged_map_procrustes/fig_s63_individual_maps_merged_map_procrustes.R")

message("======== Generating figure S64 (merged map vs individual maps map distances scatter plots) ========")
source("som_figures/fig_s64_map_distances_scatter/fig_s64_map_distances_scatter.R")
