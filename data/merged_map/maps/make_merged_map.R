
rm(list = ls())

library(Racmacs)

source('code/data_generation/load_maps_for_comparison.R')
source('code/metadata/common.R')

maps <- load_maps_for_comparison('', sr_group_standardisation = 'data/metadata/sr_group_standardisation_WT_conv_as_D614G_as614D.csv')

# Merged map with duplicated antigens only
# Make merged map from all titers
mergedD <- mergeMaps(maps)
# Subset map to antigens occurring in more than one dataset
mergedD <- removeAntigens(mergedD, agNames(mergedD)[!(agNames(mergedD) %in% duplicated_ags)])
# Optimize map
dilutionStepsize(mergedD) <- 0
mergedD <- optimizeMap(mergedD, 2, 5000)

# Warning message:
#   In optimizeMap(mergedAll, 2, 5000) :
#   There is some variation (4.33 AU for one point) in the top runs, this may be an indication that more optimization runs could help achieve a better optimum. If this still fails to help see ?unstableMaps for further possible causes.
# Realing the merged map to the duke map
mergedD <- realignMap(mergedD, maps[[1]])
# Do some styling
ptDrawingOrder(mergedD) <- rev(ptDrawingOrder(mergedD))
srOutlineWidth(mergedD) <- 1
mergedD <- rotateMap(mergedD, -11)

# Check the warning message
t <- procrustesMap(mergedD, mergedD, optimization_number = 1, comparison_optimization_number = 2, sera = F)
plot(t) # stable
t <- procrustesMap(mergedD, mergedD, optimization_number = 1, comparison_optimization_number = 3, sera = F)
plot(t) # C.36.3 moves
t <- procrustesMap(mergedD, mergedD, optimization_number = 1, comparison_optimization_number = 4, sera = F)
plot(t) # C.36.3 moves
t <- procrustesMap(mergedD, mergedD, optimization_number = 1, comparison_optimization_number = 5, sera = F)
plot(t) # C.36.3 moves

save.acmap(mergedD, 'data/merged_map/maps/merged_duplicated_antigens_only.ace', compress = T)


# Make 3D merged map with duplicated antigens only
mergedD3D <- optimizeMap(mergedD, 3, 2000)

# Align map to Duke map
mergedD3D <- realignMap(mergedD3D, maps[[1]])
# Do some styling
ptDrawingOrder(mergedD3D) <- rev(ptDrawingOrder(mergedD3D))
srOutlineWidth(mergedD3D) <- 1

save.acmap(mergedD3D, 'data/merged_map/maps/merged_duplicated_antigens_only_3D.ace', compress = T)

# Make merged map leaving each individual map out in turn, for cross validation
for (n in 1:length(maps)) {

  # Retain the map that is currently removed
  removed <- maps[[n]]

  print(mapName(removed))

  # Remove the current map
  maps_subset <- maps[-n]

  # Make merged map from subset titers
  mergedSubset <- mergeMaps(maps_subset)
  
  # Remove non-duplicated antigens
  mergedSubset <- removeAntigens(mergedSubset, agNames(mergedSubset)[!(agNames(mergedSubset) %in% duplicated_ags)])
  dilutionStepsize(mergedSubset) <- 0
  mergedSubset <- suppressMessages(optimizeMap(mergedSubset, 2, 2000, verbose = FALSE))

  mergedSubset <- realignMap(mergedSubset, maps[[1]])
  ptDrawingOrder(mergedSubset) <- rev(ptDrawingOrder(mergedSubset))

  save.acmap(mergedSubset, paste0('data/merged_map/maps/merged_duplicated_antigens_without_', mapName(removed), '.ace'), compress = T)
}


# Make merged map leaving out the CPE assay maps, for cross validation
maps_subset <- maps[-c(2, 5, 6)]
mergedSubset <- mergeMaps(maps_subset)
mergedSubset <- removeAntigens(mergedSubset, agNames(mergedSubset)[!(agNames(mergedSubset) %in% duplicated_ags)])
dilutionStepsize(mergedSubset) <- 0
mergedSubset <- optimizeMap(mergedSubset, 2, 2000, verbose = FALSE)

# Warning message:
#   In optimizeMap(mergedSubset, 2, 2000, verbose = FALSE) :
#   There is some variation (7.55 AU for one point) in the top runs, this may be an indication that more optimization runs could help achieve a better optimum. If this still fails to help see ?unstableMaps for further possible causes.

# Do some styling 
mergedSubset <- realignMap(mergedSubset, maps[[1]])
ptDrawingOrder(mergedSubset) <- rev(ptDrawingOrder(mergedSubset))
srOutlineWidth(mergedSubset) <- 1
mergedSubset <- rotateMap(mergedSubset, -7)

# Check the warning message
t <- procrustesMap(mergedSubset, mergedSubset, optimization_number = 1, comparison_optimization_number = 2, sera = F)
plot(t) # stable
t <- procrustesMap(mergedSubset, mergedSubset, optimization_number = 1, comparison_optimization_number = 3, sera = F)
plot(t) # BA.5 moves ~0.5 units
t <- procrustesMap(mergedSubset, mergedSubset, optimization_number = 1, comparison_optimization_number = 4, sera = F)
plot(t) # stable
t <- procrustesMap(mergedSubset, mergedSubset, optimization_number = 1, comparison_optimization_number = 5, sera = F)
plot(t) # stable

save.acmap(mergedSubset, paste0('data/merged_map/maps/merged_duplicated_antigens_without_cpe.ace'), compress = T)

