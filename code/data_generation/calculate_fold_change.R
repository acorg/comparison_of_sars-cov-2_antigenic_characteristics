#' Calculate a table of fold change for each map, serum group and antigen.
#' @param map_info The titers.
#' @param sr_groups The serum groups for which fold change should be calculated.
#' @param homologous_info The homologous antigens for each map and serum group
#' @export
calculate_fold_change <- function(map_info, sr_groups, homologous_info) {
  foldchange_table <- tibble(
    ag_name = character(0),
    sr_group = character(0),
    map = character(0),
    mean_diff = numeric(0),
    mean_diff_lower = numeric(0),
    mean_diff_upper = numeric(0)
  )
  
  # Compute fold drops for map
  for (sg in sr_groups) {
    for (map_name in unique(map_info$map)) {
      homologous_ag <- unname(homologous_ags[paste(sg, map_name)])
      homologous_titers <- subset(map_info, sr_group == sg & ag_name == homologous_ag & map == map_name)$titer
        
      for (ag in unique(subset(map_info, sr_group == sg & map == map_name)$ag_name)) {
        print(paste(sg, map_name, ag, homologous_ag))
        if (ag == homologous_ag) {
          if (length(setdiff(unique(homologous_titers), c('.'))) == 0) {
            # The homologous antigen didn't get titrated, add NA
            print('is homologous and doesnt exist difference NA')
            foldchange_table <- bind_rows(
              foldchange_table,
              tibble(
                ag_name = ag,
                sr_group = sg,
                map = map_name,
                mean_diff = NA,
                mean_diff_lower = NA,
                mean_diff_upper = NA
              )
            )
          } else {
            # The homologous ag did get titrated, the difference will be 0, therefore, add mean_diff = 0.
            print('is homologous, exists, difference 0')
            foldchange_table <- bind_rows(
              foldchange_table,
              tibble(
                ag_name = ag,
                sr_group = sg,
                map = map_name,
                mean_diff = 0,
                mean_diff_lower = NA,
                mean_diff_upper = NA
              )
            )
          }
        }
        if (ag != homologous_ag) {
          titers1 <- subset(map_info, sr_group == sg & ag_name == ag & map == map_name)$titer

          titers <- tibble(titers1, homologous_titers)
          titers <- subset(titers, !(titers1 %in% c('.', '*')) & !(homologous_titers %in% c('.', '*')))
          print(titers$homologous_titers)
          print(titers$titers1)

          if (length(titers$homologous_titers) >= 1) {
            # Both the homologous and the comparison antigen got titrated
            print('Both the homologous and the comparison antigen got titrated, calculating')
            fold_drop <- titertools::log2diff(
              titers2 = titers$titers1,
              titers1 = titers$homologous_titers,
              dilution_stepsize = 0,
              ci_method = 'HDI',
              mu_prior_mu = 0,
              mu_prior_sigma = 100,
              sigma_prior_alpha = 2,
              sigma_prior_beta = 0.75
            )

            foldchange_table <- bind_rows(
              foldchange_table,
              tibble(
                ag_name = ag,
                sr_group = sg,
                map = map_name,
                mean_diff = as.numeric(fold_drop[1]),
                mean_diff_lower = as.numeric(fold_drop[3]),
                mean_diff_upper = as.numeric(fold_drop[5])
              )
            )
          } else {
            # The comparison antigen did not get titrated
            print('The comparison antigen did not get titrated, adding NA')
            foldchange_table <- bind_rows(
              foldchange_table,
              tibble(
                ag_name = ag,
                sr_group = sg,
                map = map_name,
                mean_diff = NA,
                mean_diff_lower = NA,
                mean_diff_upper = NA
              )
            )
          }
        }
      }
    }
  }
  foldchange_table
}


#' Calculate a table of fold change for each map, serum group and antigen, with
#' an x-axis offset for plotting.
#' @param map_info The titers.
#' @param sr_groups The serum groups for which fold change should be calculated.
#' @param homologous_info The homologous antigens for each map and serum group
#' @param datasetorder Whether the datasets should be ordered by animal model or
#' by assay.
#' @export
calculate_xloc_foldchange_table <- function(
  foldchange_table, mean_fold_change, sr_groups,
  datasetorder = 'animal') {

  # Make table to look up the x-position of each antigen
  x_pos_lookup <- tibble(
    ag_name = character(0),
    sr_group = character(0),
    gmt = numeric(0),
    gmt_titer = numeric(0),
    xx = numeric(0),
    toSearch = character(0)
  )


  for(sg in sr_groups) {
    y <- subset(mean_fold_change, sr_group == sg)

    y[order(y$mean_, decreasing = T),] %>%
      ungroup() %>%
      mutate(
        xx = seq(0, 0.95, by=0.95 / 22),
        toSearch = paste(sr_group, ag_name)
      ) -> yy

    x_pos_lookup <- bind_rows(x_pos_lookup, unique(yy))

  }
  
  foldchange_table %>%
    ungroup()  %>%
    mutate(
      map = factor(map, levels = animal_dataset_order)
    ) %>%
    mutate(
      ag_name = as.character(ag_name)
    ) %>%
    complete(
      map, sr_group, ag_name
    ) %>%
    arrange(
      sr_group, map, desc(mean_diff)
    ) %>%
    filter(sr_group %in% sr_groups) %>%
    filter(ag_name %in% duplicated_ags) %>%
    mutate(
      sr_group = factor(sr_group, levels = c(sr_groups))
    )  -> x

  # Add the correct x location to the dataframe
  if (datasetorder == 'animal') {
    x_loc_fct <- function(sg, ag, m) {

      xpos <- subset(x_pos_lookup, toSearch == paste(sg, ag))$xx
      
      if (m == 'duke') {
        xpos <- xpos + 0
      } else if (m == 'emory')
      {
        xpos <- xpos + 1
      } else if (m == 'fda')
      {
        xpos <- xpos + 2
      } else if (m == 'innsbruck')
      {
        xpos <- xpos + 3
      } else if (m == 'oxford')
      {
        xpos <- xpos + 4
      } else if (m == 'mt_sinai_human')
      {
        xpos <- xpos + 5
      } else if (m == 'amc')
      {
        xpos <- xpos + 6
      } else if (m == 'geneva')
      {
        xpos <- xpos + 7
      } else if (m == 'emc_prnt')
      {
        xpos <- xpos + 9
      } else if (m == 'charite')
      {
        xpos <- xpos + 8
      } else if (m == 'emc_vero')
      {
        xpos <- xpos + 10
      } else if (m == 'emc_calu')
      {
        xpos <- xpos + 11
      } else if (m == 'galveston')
      {
        xpos <- xpos + 12
      } else if (m == 'madison_frnt')
      {
        xpos <- xpos + 13
      } else if (m == 'madison_pooled')
      {
        xpos <- xpos + 14
      } else if (m == 'madison_unpooled')
      {
        xpos <- xpos + 15
      } else if (m == 'st_louis')
      {
        xpos <- xpos + 16
      } else if (m == 'maryland')
      {
        xpos <- xpos + 17
      }
      xpos
    }
    } else {
      x_loc_fct <- function(sg, ag, m) {
        
        xpos <- subset(x_pos_lookup, toSearch == paste(sg, ag))$xx
        
        if (m == 'emory') {
          xpos <- xpos + 0
        } else if (m == 'innsbruck')
        {
          xpos <- xpos + 1
        } else if (m == 'oxford')
        {
          xpos <- xpos + 2
        } else if (m == 'galveston')
        {
          xpos <- xpos + 3
        } else if (m == 'st_louis')
        {
          xpos <- xpos + 5
        } else if (m == 'madison_frnt')
        {
          xpos <- xpos + 4
        } else if (m == 'duke')
        {
          xpos <- xpos + 6
        } else if (m == 'fda')
        {
          xpos <- xpos + 7
        } else if (m == 'amc')
        {
          xpos <- xpos + 8
        } else if (m == 'emc_vero')
        {
          xpos <- xpos + 9
        } else if (m == 'emc_calu')
        {
          xpos <- xpos + 10
        } else if (m == 'geneva')
        {
          xpos <- xpos + 11
        } else if (m == 'charite')
        {
          xpos <- xpos + 12
        } else if (m == 'emc_prnt')
        {
          xpos <- xpos + 13
        } else if (m == 'mt_sinai_human')
        {
          xpos <- xpos + 14
        } else if (m == 'madison_pooled')
        {
          xpos <- xpos + 15
        } else if (m == 'madison_unpooled')
        {
          xpos <- xpos + 16
        } else if (m == 'maryland')
        {
          xpos <- xpos + 17
        }
        xpos
    }
    
  }
  
  x$x_loc <- mapply(x_loc_fct, x$sr_group, x$ag_name, x$map)
  

  # Set outlier antigens (BA.1 to 0), all titers are ND.
  x$mean_diff[x$ag_name == 'BA.1' & x$map == 'emory' & x$sr_group == 'D614G convalescent'] <- NA
  x$mean_diff_lower[x$ag_name == 'BA.1' & x$map == 'emory' & x$sr_group == 'D614G convalescent'] <- NA
  x$mean_diff_upper[x$ag_name == 'BA.1' & x$map == 'emory' & x$sr_group == 'D614G convalescent'] <- NA

  x$mean_diff[x$ag_name == 'BA.1' & x$map == 'innsbruck' & x$sr_group == 'B.1.1.7 convalescent'] <- NA
  x$mean_diff_lower[x$ag_name == 'BA.1' & x$map == 'innsbruck' & x$sr_group == 'B.1.1.7 convalescent'] <- NA
  x$mean_diff_upper[x$ag_name == 'BA.1' & x$map == 'innsbruck' & x$sr_group == 'B.1.1.7 convalescent'] <- NA

  x

}

#' For the sampled data, make a dataframe with the correct x-location.
#' @param data_ag: The information about the estimated fold change per antigen and serum group
#' @param data_eff: The information about the estimated slope effects.
#' @param sr_groups: The serum groups which should be included.
#' @param datasetorder: The way the datasets should be ordered. Either 'animal' or 'assay'.
calculate_xloc_sampling <- function(data_ag, data_eff, titrated, sr_groups, datasetorder = 'animal') {
  data_loc <- tibble(
    variable = character(0),
    ci_low = numeric(0),
    ci_high = numeric(0),
    mean_ = numeric(0),
    ag_name = character(0),
    sr_group = character(0),
    titrated = character(0),
    foldchange = numeric(0),
    map = character(0),
    map_x_loc = numeric(0),
    group = factor(0)
  )

  # Add a column with the index of each antigen by map and serum group.
  for (map_ in unique(data_eff$map)) {
    data_ag$foldchange <- unname(data_ag$mean_ * subset(data_eff, map == map_)$mean_)
    data_ag$map <- map_
    for (srg in sr_groups) {
      d <- filter(data_ag, sr_group == srg)
      data_ag_or <- d[order(d$foldchange,decreasing=TRUE),]

      if (paste(map_, srg) %in% titrated) {
        data_ag_or$map_x_loc <- seq(0, 0.95, by=0.95 / 22)
      } else {
        data_ag_or$map_x_loc <- rep(NA, 23)
      }
      data_ag_or$group <- factor(c(c(1), rep(2, 21), c(1)))

      # Add a row with an empty antigen, so that the line breaks correctly
      data_ag_or <- bind_rows(data_ag_or, tibble(
        variable = '',
        ci_low = NA,
        ci_high = NA,
        mean_ = NA,
        ag_name = 'x',
        sr_group = srg,
        titrated = NA,
        foldchange = NA,
        map = map_,
        map_x_loc = 0.96,
        group = factor(1),
      ))

      data_loc <- bind_rows(data_loc, data_ag_or)
    }
  }

  if (datasetorder == 'animal') {
    x_loc_fct <- function(sg, ag, m) {

      xpos <- subset(data_loc, sr_group == sg & ag_name == ag & map == m)$map_x_loc

      if (m == 'duke') {
        xpos <- xpos + 0
      } else if (m == 'emory')
      {
        xpos <- xpos + 1
      } else if (m == 'fda')
      {
        xpos <- xpos + 2
      } else if (m == 'innsbruck')
      {
        xpos <- xpos + 3
      } else if (m == 'oxford')
      {
        xpos <- xpos + 4
      } else if (m == 'mt_sinai_human')
      {
        xpos <- xpos + 5
      } else if (m == 'amc')
      {
        xpos <- xpos + 6
      } else if (m == 'geneva')
      {
        xpos <- xpos + 7
      } else if (m == 'emc_prnt')
      {
        xpos <- xpos + 9
      } else if (m == 'charite')
      {
        xpos <- xpos + 8
      } else if (m == 'emc_vero')
      {
        xpos <- xpos + 10
      } else if (m == 'emc_calu')
      {
        xpos <- xpos + 11
      } else if (m == 'galveston')
      {
        xpos <- xpos + 12
      } else if (m == 'madison_frnt')
      {
        xpos <- xpos + 13
      } else if (m == 'madison_pooled')
      {
        xpos <- xpos + 14
      } else if (m == 'madison_unpooled')
      {
        xpos <- xpos + 15
      } else if (m == 'st_louis')
      {
        xpos <- xpos + 16
      } else if (m == 'maryland')
      {
        xpos <- xpos + 17
      }
      xpos
    }
  } else {
    x_loc_fct <- function(sg, ag, m) {

      xpos <- subset(data_loc, sr_group == sg & ag_name == ag & map == m)$map_x_loc

      if (m == 'emory') {
        xpos <- xpos + 0
      } else if (m == 'innsbruck')
      {
        xpos <- xpos + 1
      } else if (m == 'oxford')
      {
        xpos <- xpos + 2
      } else if (m == 'galveston')
      {
        xpos <- xpos + 3
      } else if (m == 'st_louis')
      {
        xpos <- xpos + 5
      } else if (m == 'madison_frnt')
      {
        xpos <- xpos + 4
      } else if (m == 'duke')
      {
        xpos <- xpos + 6
      } else if (m == 'fda')
      {
        xpos <- xpos + 7
      } else if (m == 'amc')
      {
        xpos <- xpos + 8
      } else if (m == 'emc_vero')
      {
        xpos <- xpos + 9
      } else if (m == 'emc_calu')
      {
        xpos <- xpos + 10
      } else if (m == 'geneva')
      {
        xpos <- xpos + 11
      } else if (m == 'charite')
      {
        xpos <- xpos + 12
      } else if (m == 'emc_prnt')
      {
        xpos <- xpos + 13
      } else if (m == 'mt_sinai_human')
      {
        xpos <- xpos + 14
      } else if (m == 'madison_pooled')
      {
        xpos <- xpos + 15
      } else if (m == 'madison_unpooled')
      {
        xpos <- xpos + 16
      } else if (m == 'maryland')
      {
        xpos <- xpos + 17
      }
      xpos
    }

  }

  data_loc$x_loc <- mapply(x_loc_fct, data_loc$sr_group, data_loc$ag_name, data_loc$map)

  data_loc
}
