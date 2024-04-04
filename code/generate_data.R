
### GENERATE DATA

## Generating individual maps

message("======== Generating duke map ========")
message("For details, check https://github.com/acorg/mapping_SARS-CoV-2_antigenic_relationships_and_serological_responses/")

message("======== Generating emory map ========")
source("data/individual_datasets/emory/maps/making_the_maps.R")
source("data/individual_datasets/emory/outliers/determine_and_remove_outliers.R")
source("data/individual_datasets/emory/individual_effects/calculate_individual_effects.R")

message("======== Generating fda map ========")
source("data/individual_datasets/fda/maps/making_the_maps.R")
source("data/individual_datasets/fda/outliers/determine_and_remove_outliers.R")
source("data/individual_datasets/fda/individual_effects/calculate_individual_effects.R")

message("======== Generating innsbruck map ========")
source("data/individual_datasets/innsbruck/maps/making_the_maps.R")
source("data/individual_datasets/innsbruck/outliers/determine_and_remove_outliers.R")
source("data/individual_datasets/innsbruck/individual_effects/calculate_individual_effects.R")

message("======== Generating oxford map ========")
source("data/individual_datasets/oxford/maps/making_the_maps.R")
source("data/individual_datasets/oxford/outliers/determine_and_remove_outliers.R")
source("data/individual_datasets/oxford/individual_effects/calculate_individual_effects.R")

message("======== Generating mt_sinai map ========")
source("data/individual_datasets/mt_sinai/maps/making_the_maps.R")
source("data/individual_datasets/mt_sinai/outliers/determine_and_remove_outliers.R")
source("data/individual_datasets/mt_sinai/individual_effects/calculate_individual_effects.R")

message("======== Generating amc map ========")
source("data/individual_datasets/amc/maps/making_the_maps.R")
source("data/individual_datasets/amc/outliers/determine_and_remove_outliers.R")
source("data/individual_datasets/amc/individual_effects/calculate_individual_effects.R")

message("======== Generating geneva map ========")
source("data/individual_datasets/geneva/maps/making_the_maps.R")
source("data/individual_datasets/geneva/outliers/determine_and_remove_outliers.R")
source("data/individual_datasets/geneva/individual_effects/calculate_individual_effects.R")

message("======== Generating madison_pooled map ========")
source("data/individual_datasets/madison_pooled/maps/making_the_maps.R")

message("======== Generating charite map ========")
source("data/individual_datasets/charite/maps/making_the_maps.R")

message("======== Generating emc_prnt, emc_vero, emc_calu map ========")
source("data/individual_datasets/emc/maps/making_the_maps.R")

message("======== Generating madison_unpooled map ========")
source("data/individual_datasets/madison_unpooled/maps/making_the_maps.R")

message("======== Generating galveston map ========")
source("data/individual_datasets/galveston/maps/making_the_maps.R")

message("======== Generating madison_frnt map ========")
source("data/individual_datasets/madison_frnt/maps/making_the_maps.R")

message("======== Generating maryland_pooled map ========")
source("data/individual_datasets/maryland/maps/making_the_maps.R")

message("======== Generating wustl map ========")
source("data/individual_datasets/wustl/maps/making_the_maps.R")


## Generating titer data

# Generating titer mangitude related data
message("======== Generating titer magnitude effect data ========")
source("data/titer_analyses/titer_magnitude/dataset_magnitude_effect/estimate_dataset_magnitude_effect.R")
source("data/titer_analyses/titer_magnitude/dataset_magnitude_effect/dataset_magnitude_effect_sampling.R")

message("======== Titer magnitude effect by species and assay ========")
source("data/titer_analyses/titer_magnitude/dataset_magnitude_effect_animal_model_assay_differences/estimate_dataset_magnitude_effect_animal_model_assay_differences_sampling.R")

message("======== Titer magnitude effect by species and assay without non-NT50 datasets ========")
source("data/titer_analyses/titer_magnitude/dataset_magnitude_effect_animal_model_assay_differences_nt50/estimate_dataset_magnitude_effect_animal_model_assay_differences_nt50_sampling.R")

message("======== Titer magnitude effect by titer unit ========")
source("data/titer_analyses/titer_magnitude/dataset_magnitude_effect_titer_units/estimate_dataset_magnitude_effect_titer_magnitude_titer_units_sampling.R")

# Generating titer variability related data
message("======== Generating titer variability data ========")
source("data/titer_analyses/titer_variability/titer_variability_effect/summarise_titer_variability.R")

# Generating fold change data
message("======== Generating fold change data ========")
source("data/titer_analyses/foldchange/fold_change_calculation/fold_change_calculation.R")
source("data/titer_analyses/foldchange/fold_change_kendall/fold_change_kendall.R")

message("======== Generating fold change slope data ========")
source("data/titer_analyses/foldchange/slope_calculation/slope_calculation.R")

message("======== Generating fold change slope data by animal and assay ========")
source("data/titer_analyses/foldchange/slope_calculation_by_animal_assay/slope_calculation_by_animal_assay.R")

# Generating immunodominance data
message("======== Generating immunodominance data ========")
source("data/titer_analyses/immunodominance/immunodominance_data/make_immunodominance_data.R")


## Generating the merged map

message("======== Generating merged maps ========")
# This generates the following maps:
#   The merged map using all individual maps
#   The merged map using all individual maps but only including antigens ocurring in at least two individual maps
#   The merged map optimised in 3 dimensions using all individual maps but only including antigens ocurring in at least two individual maps
#   18 merged maps, leaving each individual map out in turn
#   A merged map without individual maps made with CPE assay
source("data/merged_map/maps/make_merged_map.R")
