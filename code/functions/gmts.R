
#' Summarise the GMTs
summarise_gmts <- function(mapdata) {

  mapdata %>%
    summarise(
      gmt_stat = list(
        titertools::gmt(
          titers = titer,
          ci_method = 'HDI',
          sigma_prior_alpha = 2,
          sigma_prior_beta = 0.75,
          mu_prior_mu = 0,
          mu_prior_sigma = 100,
          dilution_stepsize = unique(dilution_stepsize)
        )
      ),
      .groups = 'drop'
    ) %>%
    mutate(
      gmt = vapply(gmt_stat, \(x) x['mean', 'estimate'], numeric(1)),
      gmt_upper = vapply(gmt_stat, \(x) x['mean', 'lower'], numeric(1)),
      gmt_lower = vapply(gmt_stat, \(x) x['mean', 'upper'], numeric(1))
    )
  
}


#' Summarise the GMTs for each serum group and variant in a map
#' @param map An acmap anitgenic map
#' @param ci_method The method to use to infer the confidence interval. Must be one
#' of 'quap', 'ETI', 'HDI', 'BCI', or 'SI'.
#' @export
srGroupGMTs <- function(map, ci_method = 'HDI') {

  # Fetch titer table
  titer_table <- adjustedTiterTable(map)

  # Calculate gmts for each group and antigen
  sr_group_gmts <- lapply(
    unique(srGroups(map)),
    \(sr_group) {
      apply(
        titer_table[ , srGroups(map) == sr_group], 1, \(titers) {
          titertools::gmt(
            titers = titers,
            ci_method = ci_method,
            sigma_prior_alpha = 2,
            sigma_prior_beta = 0.75,
            mu_prior_mu = 0,
            mu_prior_sigma = 100,
            dilution_stepsize = dilutionStepsize(map)
          )['mean', 'estimate']
        }
      )
    }
  )

  # Convert to a matrix and return the values
  sr_group_gmts <- do.call(cbind, sr_group_gmts)
  colnames(sr_group_gmts) <- as.character(unique(srGroups(map)))
  sr_group_gmts <- sr_group_gmts[ ,match(levels(srGroups(map)), colnames(sr_group_gmts))]
  sr_group_gmts

}
