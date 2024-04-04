# '# Calculate mean fold change difference for eacg datasets
#' @param map_info The titers.
#' @param sr_groups The serum groups for which fold change should be calculated.
summarise_map_folddrop <- function(map_info, sr_groups) {
  curated_data_aggregated_xxx_animals <- tibble(
    sr_group = character(0),
    titer_diff = numeric(0),
    titer_diff_upper = numeric(0),
    titer_diff_lower = numeric(0),
    position = character(0),
    map = character(0)
  )

  for (map_ in unique(map_info$map)) {
    for (sg in sr_groups) {
      titers1 <- c()
      titers2 <- c()
      for (substitution_pair in subst_pairs_484) {
        subst1 <- substitution_pair[1]
        subst2 <- substitution_pair[2]

        titers1map <- subset(map_info, map == map_ & sr_group == sg & ag_name == subst1)$titer
        titers2map <- subset(map_info, map == map_ & sr_group == sg & ag_name == subst2)$titer

        titers1 <- c(titers1, titers1map)
        titers2 <- c(titers2, titers2map)
      }


      titersMap <- tibble(titers1, titers2)
      titersMap <- subset(titersMap, !(titers1 %in% c('.', '*')) & !(titers2 %in% c('.', '*')))

      if (length(titersMap$titers1) >= 1 & length(titersMap$titers2) >= 1) {
        result <- titertools::log2diff(
          titers1 = titersMap$titers1,
          titers2 = titersMap$titers2,
          dilution_stepsize = 0,
          ci_method = 'HDI',
          sigma_prior_alpha = 2,
          sigma_prior_beta = 0.75
        )

        curated_data_aggregated_xxx_animals <- rbind(curated_data_aggregated_xxx_animals, tibble(
          sr_group = sg,
          titer_diff = result[1],
          titer_diff_upper = result[5],
          titer_diff_lower = result[3],
          position = '484',
          map = map_
        )
        )

      } else {
        curated_data_aggregated_xxx_animals <- rbind(curated_data_aggregated_xxx_animals, tibble(
          sr_group = sg,
          titer_diff = NA,
          titer_diff_upper = NA,
          titer_diff_lower = NA,
          position = '484',
          map = map_
        )
        )
      }
    }
  }

  for (map_ in unique(map_info$map)) {
    for (sg in sr_groups) {
      titers1 <- c()
      titers2 <- c()
      for (substitution_pair in subst_pairs_501) {
        subst1 <- substitution_pair[1]
        subst2 <- substitution_pair[2]

        titers1map <- subset(map_info, map == map_ & sr_group == sg & ag_name == subst1)$titer
        titers2map <- subset(map_info, map == map_ & sr_group == sg & ag_name == subst2)$titer

        titers1 <- c(titers1, titers1map)
        titers2 <- c(titers2, titers2map)
      }


      titersMap <- tibble(titers1, titers2)
      titersMap <- subset(titersMap, !(titers1 %in% c('.', '*')) & !(titers2 %in% c('.', '*')))

      if (length(titersMap$titers1) >= 1 & length(titersMap$titers2) >= 1) {
        result <- titertools::log2diff(
          titers1 = titersMap$titers1,
          titers2 = titersMap$titers2,
          dilution_stepsize = 0,
          ci_method = 'HDI',
          sigma_prior_alpha = 2,
          sigma_prior_beta = 0.75
        )

        curated_data_aggregated_xxx_animals <- rbind(curated_data_aggregated_xxx_animals, tibble(
          sr_group = sg,
          titer_diff = result[1],
          titer_diff_upper = result[5],
          titer_diff_lower = result[3],
          position = '501',
          map = map_
        )
        )

      } else {
        curated_data_aggregated_xxx_animals <- rbind(curated_data_aggregated_xxx_animals, tibble(
          sr_group = sg,
          titer_diff = NA,
          titer_diff_upper = NA,
          titer_diff_lower = NA,
          position = '501',
          map = map_
        )
        )
      }
    }
  }

  curated_data_aggregated_xxx_animals %>%
    mutate(
      sr_group = factor(sr_group, levels = sr_groups)
    ) -> curated_data_aggregated_xxx_animals

  curated_data_aggregated_xxx_animals
}
