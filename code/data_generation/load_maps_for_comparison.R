
#' @export
standardise_sr_groups <- function(map, standardisation_info) {
  sr_groups <- standardisation_info$Standardisation[match(srGroups(map), standardisation_info$Serum)]
  srGroups(map) <- factor(sr_groups, levels = unique(sr_groups))
  map
}

#' A function that loads a set of desired maps and performs standardisation
#' on antigen and serum names.
#' @param base_path The path to the directory where you store your SARS2 repos.
#' @param sr_group_standardisation A file with the serum group standardisation.
#' @export
load_maps_for_comparison <- function(
    base_path,
    sr_group_standardisation = NULL
    ){

    # Duke
    duke_path <- paste0(base_path, 'data/individual_datasets/duke/maps/map_no_outliers.ace')
    duke <- read.acmap(duke_path)
    mapName(duke) <- 'duke'

    agNames(duke)[agNames(duke) == 'BA.4/BA.5'] <- 'BA.5'
    agSize(duke)[agNames(duke) == 'BA.5'] <- 18
    agSize(duke)[agNames(duke) == 'BA.2'] <- 18
    agFill(duke)[agNames(duke) == 'BA.5'] <- '#F08DA5'
    agFill(duke)[agNames(duke) == 'BA.2'] <- '#d10fa2'
    agFill(duke)[agNames(duke) == 'BA.2.12.1'] <- '#d10fa2'

    srGroups(duke) <- forcats::fct_drop(
        forcats::fct_collapse(
            srGroups(duke),
            'D614G convalescent' = c('D614G'),
            'B.1.1.7 convalescent' = c('B.1.1.7'),
            'P.1 convalescent' = c('P.1'),
            'B.1.351 convalescent' = c('B.1.351'),
            'B.1.526+E484K convalescent' = c('B.1.526+E484K'),
            'B.1.617.2 convalescent' = c('B.1.617.2'),
            'B.1.637 convalescent' = c('B.1.637'),
            'C.37 convalescent' = c('C.37'),
            'BA.1 convalescent' = c('BA.1'),
            'BA.2 convalescent' = c('BA.2'),
            'mRNA-1273' = c('2x mRNA-1273'),
            '3x mRNA-1273 BD01'   = c('3x mRNA-1273 BD01'  ),
            '3x mRNA-1273 BD29' = c('3x mRNA-1273 BD29'),
            '3x mRNA-1273 (6 month)' = c('3x mRNA-1273 (6 month)'),
            'mRNA-1273.351' = c('2x mRNA-1273.351')
        )
    )

    # Change point drawing order so important ags are at the front
    bring_to_front <- c('B.1.617.2', 'BA.1', 'D614G', 'B.1.351')
    indices <- which(agNames(duke) %in% bring_to_front)
    ptDrawingOrder(duke) <- c(
      ptDrawingOrder(duke)[!ptDrawingOrder(duke) %in% indices],
      indices)

    duke <- rotateMap(duke, -30)
    duke <- translateMap(duke, c(0.7, -0.73))

    # Maryland
    maryland_path <- paste(base_path,
        'data/individual_datasets/maryland/maps/map.ace', sep='')
    maryland <- read.acmap(maryland_path)

    maryland <- rotateMap(maryland, -65.5)
    maryland <- translateMap(maryland, c(2.4, -1.1))

    # Galveston 28 days
    galveston_path <- paste(base_path,
        'data/individual_datasets/galveston/maps/map.ace', sep='')
    galveston <- read.acmap(galveston_path)

    # Change point drawing order so important ags are at the front
    bring_to_front <- c('614D', 'B.1.617.2', 'B.1.351')
    indices <- which(agNames(galveston) %in% bring_to_front)
    ptDrawingOrder(galveston) <- c(
        ptDrawingOrder(galveston)[!ptDrawingOrder(galveston) %in% indices],
        indices)

    galveston <- translateMap(galveston, c(0.1, 0))


    # Emory
    emory_path <- paste0(base_path, 'data/individual_datasets/emory/maps/map_no_outliers.ace')
    emory <- read.acmap(emory_path)
    agFill(emory)[agNames(emory) == 'B.1.630'] <- '#9dd455'

    # Change point drawing order so important ags are at the front
    bring_to_front <- c('B.1.1.7', 'B.1.351', 'P.1', 'B.1.526+E484K',
        'B.1.525', 'P.2', 'D614G', '614D')
    indices <- which(agNames(emory) %in% bring_to_front)
    ptDrawingOrder(emory) <- c(
        ptDrawingOrder(emory)[!ptDrawingOrder(emory) %in% indices], indices)

    emory <- translateMap(emory, c(0.3, 0))
    emory <- rotateMap(emory, 7)


    # Madison pooled
    mad_path <- paste(base_path,
        'data/individual_datasets/madison_pooled/maps/map.ace', sep='')
    mad_pooled <- read.acmap(mad_path)
    titerTable(mad_pooled) <- titerTable(mad_pooled) # Hack to get around Racmacs error merging tables with multiple titer layers

    # Change point drawing order so important ags are at the front
    ptDrawingOrder(mad_pooled) <- rev(ptDrawingOrder(mad_pooled))

    bring_to_front <- c('D614G')
    indices <- which(agNames(mad_pooled) %in% bring_to_front)
    ptDrawingOrder(mad_pooled) <- c(
        ptDrawingOrder(mad_pooled)[!ptDrawingOrder(mad_pooled) %in% indices],
        indices)

    mad_pooled <- reflectMap(mad_pooled, axis = 'x')
    mad_pooled <- rotateMap(mad_pooled, -120)
    mad_pooled <- translateMap(mad_pooled, c(1.8, -1))

    # Madison unpooled
    mad_unpooled_path <- paste(base_path,
        'data/individual_datasets/madison_unpooled/maps/map.ace',
        sep='')

    mad_unpooled <- read.acmap(mad_unpooled_path)
    titerTable(mad_unpooled) <- titerTable(mad_unpooled)

    ptDrawingOrder(mad_unpooled) <- rev(ptDrawingOrder(mad_unpooled))

    mapName(mad_unpooled) <- 'madison_unpooled'

    mad_unpooled <- rotateMap(mad_unpooled, -10)
    mad_unpooled <- translateMap(mad_unpooled, c(2.2, -1.1))


    # Madison FRNT
    mad_frnt_path <- paste(base_path,
        'data/individual_datasets/madison_frnt/maps/map.ace',
        sep='')

    mad_frnt <- read.acmap(mad_frnt_path)

    mad_frnt <- rotateMap(mad_frnt, -19)
    mad_frnt <- translateMap(mad_frnt, c(1, -1.8))
    mapName(mad_frnt) <- 'madison_frnt'

    # Washington University
    wustl_path <- paste(base_path,
        'data/individual_datasets/wustl/maps/map.ace', sep='')
    wustl <- read.acmap(wustl_path)

    # Change point drawing order so important ags are at the front
    bring_to_front <- c('B.1.351')
    indices <- which(agNames(wustl) %in% bring_to_front)
    ptDrawingOrder(wustl) <- c(
        ptDrawingOrder(wustl)[!ptDrawingOrder(wustl) %in% indices],
        indices)

    mapName(wustl) <- 'st_louis'

    wustl <- rotateMap(wustl, -34)
    wustl <- translateMap(wustl, c(0.75, -1.1))

    # Oxford
    oxford_path <- paste(base_path,
        'data/individual_datasets/oxford/maps/map_no_outliers.ace', sep='')
    oxford <- read.acmap(oxford_path)

    # Change point drawing order so important ags are at the front
    bring_to_front <- c('B.1.617.2')
    indices <- which(agNames(oxford) %in% bring_to_front)
    ptDrawingOrder(oxford) <- c(
        ptDrawingOrder(oxford)[!ptDrawingOrder(oxford) %in% indices],
        indices)

    oxford <- translateMap(oxford, c(0.45, -0.25))

    # Mt Sinai human
    mt_sinai_h_path <- paste(base_path,
        'data/individual_datasets/mt_sinai/maps/map_no_outliers.ace', sep='')
    mt_sinai_h <- read.acmap(mt_sinai_h_path)

    # Change point drawing order so important ags are at the front
    bring_to_front <- c('614D', 'B.1.617.2')
    indices <- which(agNames(mt_sinai_h) %in% bring_to_front)
    ptDrawingOrder(mt_sinai_h) <- c(
        ptDrawingOrder(mt_sinai_h)[!ptDrawingOrder(mt_sinai_h) %in% indices],
        indices)

    mt_sinai_h <- translateMap(mt_sinai_h, c(0.1, -0.7))

    # EMC
    emc_path <- paste(base_path,
        'data/individual_datasets/emc/maps/emc-prnt.ace', sep='')
    emc <- read.acmap(emc_path)
    mapName(emc) <- 'emc_prnt'

    # Change point drawing order so important ags are at the front
    bring_to_front <- c('B.1.617.2')
    indices <- which(agNames(emc) %in% bring_to_front)
    ptDrawingOrder(emc) <- c(
        ptDrawingOrder(emc)[!ptDrawingOrder(emc) %in% indices],
        indices)

    emc <- rotateMap(emc, -53)
    emc <- translateMap(emc, c(0.5, -2.3))

    # Innsbruck
    innsbruck_path <- paste(base_path,
        'data/individual_datasets/innsbruck/maps/map-OmicronI+II+III-thresholded-single_exposure-P1m1-no-nds-no-outliers.ace', sep='')
    innsbruck <- read.acmap(innsbruck_path)
    titerTable(innsbruck) <- titerTable(innsbruck)
    mapName(innsbruck) <- 'innsbruck'
    # Change point drawing order so important ags are at the front
    bring_to_front <- c('B.1.351')
    indices <- which(agNames(innsbruck) %in% bring_to_front)
    ptDrawingOrder(innsbruck) <- c(
      ptDrawingOrder(innsbruck)[!ptDrawingOrder(innsbruck) %in% indices],
      indices)
    agSize(innsbruck)[agNames(innsbruck) == 'B.1.1.7+E484K'] <- 12
    srSize(innsbruck) <- 10
    agNames(innsbruck)[agNames(innsbruck) == 'P.1.1'] <- 'P.1'

    agFill(innsbruck) <- c('#393b79', '#637939', '#637939', '#7b4173', '#e7ba52', '#d18652', '#EF3737', '#d10fa2', '#F08DA5')
    srOutline(innsbruck)[srGroups(innsbruck) == 'D614G convalescent'] <- '#333333'
    srOutline(innsbruck)[srGroups(innsbruck) == 'B.1.617.2 convalescent'] <- '#d18652'
    srOutline(innsbruck)[srGroups(innsbruck) == 'B.1.1.7 convalescent'] <- '#637939'
    srOutline(innsbruck)[srGroups(innsbruck) == 'B.1.351 convalescent'] <- '#e7ba52'
    srOutline(innsbruck)[srGroups(innsbruck) == 'mRNA-1273'] <- 'grey'
    srOutline(innsbruck)[srGroups(innsbruck) == 'AstraZeneca'] <- 'grey'
    srOutline(innsbruck)[srGroups(innsbruck) == 'AstraZeneca-Pfizer'] <- 'grey'
    srOutline(innsbruck)[srGroups(innsbruck) == 'Pfizer'] <- 'grey'
    srOutline(innsbruck)[srGroups(innsbruck) == 'BA.1 convalescent'] <- '#EF3737'
    srOutline(innsbruck)[srGroups(innsbruck) == 'BA.2 convalescent'] <- '#d10fa2'

    innsbruck <- rotateMap(innsbruck, 195)
    innsbruck <- translateMap(innsbruck, c(2.2, -0.5))

    # Charite
    charite_path <- paste(base_path,
        'data/individual_datasets/charite/maps/map.ace', sep='')
    charite <- read.acmap(charite_path)
    agNames(charite) <- c('B.1.1.7', 'BA.1', 'BA.2-12', 'BA.2', 'BA.4', 'BA.5',
                          'B.1.351', 'B.1.617.2', 'B.1.621', 'D614G',
                          'B.1+E484K')
    mapName(charite) <- 'charite'
    srSize(charite) <- 10
    bring_to_front <- c('D614G', 'BA.5')
        indices <- which(agNames(charite) %in% bring_to_front)
        ptDrawingOrder(charite) <- c(
            ptDrawingOrder(charite)[!ptDrawingOrder(charite) %in% indices],
            indices)
    agFill(charite)[agNames(charite) %in% c('BA.2', 'BA.2-12')] <- '#d10fa2'
    sr_groups <- c('BA.2 convalescent', 'BA.2 convalescent', 'BA.2 convalescent', 'B.1.1.7 convalescent', 'B.1.1.7 convalescent',
                   'B.1.351 convalescent', 'B.1.351 convalescent', 'B.1.351 convalescent',
                   'B.1.617.2 convalescent', 'B.1.617.2 convalescent', 'B.1.617.2 convalescent',
                   'BA.1 convalescent', 'BA.1 convalescent', 'BA.1 convalescent',
                   'BA.2 convalescent', 'BA.2 convalescent', 'BA.2 convalescent',
                   'B.1+E484K convalescent', 'B.1+E484K convalescent', 'B.1+E484K convalescent',
                   'D614G convalescent', 'D614G convalescent', 'D614G convalescent')
    srGroups(charite) <- factor(sr_groups, levels = unique(sr_groups))
    srOutline(charite)[srGroups(charite) == 'D614G convalescent'] <- '#333333'
    srOutline(charite)[srGroups(charite) == 'BA.2 convalescent'] <- '#d10fa2'

    charite <- rotateMap(charite, -35)
    charite <- translateMap(charite, c(3.5, -1.05))

    bring_to_front <- c('BA.2')
    indices <- which(agNames(charite) %in% bring_to_front)
    ptDrawingOrder(charite) <- c(
      ptDrawingOrder(charite)[!ptDrawingOrder(charite) %in% indices],
      indices)

    # FDA
    fda_path <- paste(base_path,
        'data/individual_datasets/fda/maps/map_no_outliers.ace', sep='')
    fda <- read.acmap(fda_path)

    fda <- translateMap(fda, c(1, -0.7))
    mapName(fda) <- 'fda'
    agSize(fda)[agNames(fda) == 'BA.1.1'] <- 12
    agSize(fda)[agNames(fda) == 'BA.2.12.1'] <- 12
    agFill(fda)[agNames(fda) == 'BA.2.12.1'] <- '#d10fa2'

    bring_to_front <- c('BA.1', 'B.1.351', 'BA.5', 'BA.2')
    indices <- which(agNames(fda) %in% bring_to_front)
    ptDrawingOrder(fda) <- c(
      ptDrawingOrder(fda)[!ptDrawingOrder(fda) %in% indices],
      indices)

    fda <- rotateMap(fda, 2)
    fda <- translateMap(fda, c(0.3, 0))

    # Geneva
    geneva_path <- paste(base_path,
        'data/individual_datasets/geneva/maps/map_no_outliers.ace', sep='')
    geneva <- read.acmap(geneva_path)

    geneva <- translateMap(geneva, c(0.5, -0.25))
    mapName(geneva) <- 'geneva'

    # AMC
    amc_path <- paste(base_path,
        'data/individual_datasets/amc/maps/map_no_outliers.ace', sep='')
    amc <- read.acmap(amc_path)

    amc <- translateMap(amc, c(1, -2.1))
    amc <- rotateMap(amc, -5.5)
    mapName(amc) <- 'amc'

    # EMC Calu
    emc_calu_path <- paste(base_path,
        'data/individual_datasets/emc/maps/pseudo-calu.ace', sep='')
    emc_calu <- read.acmap(emc_calu_path)

    emc_calu <- translateMap(emc_calu, c(1.5, -1.6))
    emc_calu <- rotateMap(emc_calu, -1)
    mapName(emc_calu) <- 'emc_calu'

    # EMC Vero
    emc_vero_path <- paste(base_path,
        'data/individual_datasets/emc/maps/pseudo-vero.ace', sep='')
    emc_vero <- read.acmap(emc_vero_path)

    emc_vero <- translateMap(emc_vero, c(1.8, -1.75))
    mapName(emc_vero) <- 'emc_vero'


    if (!(is.null(sr_group_standardisation))) {
        standardisation_info <- read.csv(sr_group_standardisation, stringsAsFactors = FALSE)
        duke <- standardise_sr_groups(duke, standardisation_info)
        maryland <- standardise_sr_groups(maryland, standardisation_info)
        galveston <- standardise_sr_groups(galveston, standardisation_info)
        emory <- standardise_sr_groups(emory, standardisation_info)
        # For madison_pooled, we want to differentiate between the D614G and 614D serum
        mad_pooled <- standardise_sr_groups(mad_pooled, standardisation_info)
        mad_unpooled <- standardise_sr_groups(mad_unpooled, standardisation_info)
        wustl <- standardise_sr_groups(wustl, standardisation_info)
        oxford <- standardise_sr_groups(oxford, standardisation_info)
        mt_sinai_h <- standardise_sr_groups(mt_sinai_h, standardisation_info)
        emc <- standardise_sr_groups(emc, standardisation_info)
        innsbruck <- standardise_sr_groups(innsbruck, standardisation_info)
        charite <- standardise_sr_groups(charite, standardisation_info)
        mad_frnt <- standardise_sr_groups(mad_frnt, standardisation_info)
        fda <- standardise_sr_groups(fda, standardisation_info)
        geneva <- standardise_sr_groups(geneva, standardisation_info)
        amc <- standardise_sr_groups(amc, standardisation_info)
        emc_calu <- standardise_sr_groups(emc_calu, standardisation_info)
        emc_vero <- standardise_sr_groups(emc_vero, standardisation_info)
    }

    {
        maps <- list(
          duke = duke,
          maryland = maryland,
          galveston = galveston,
          emory = emory,
          madison_pooled = mad_pooled,
          madison_unpooled = mad_unpooled,
          st_louis = wustl,
          oxford = oxford,
          mt_sinai_human = mt_sinai_h,
          emc_prnt = emc,
          innsbruck = innsbruck,
          charite = charite,
          madison_frnt = mad_frnt,
          fda = fda,
          geneva = geneva,
          amc = amc,
          emc_calu = emc_calu,
          emc_vero = emc_vero
        )
    }
}



#' @export
groupAgs <- function(map, groupings) {

  ag_names <- agNames(map)

  for (n in seq_along(groupings)) {
    ag_names[ag_names %in% groupings[[n]]] <- names(groupings)[n]
  }

  agNames(map) <- ag_names
  map

}


#' Arrange the data for estimating the titers
#' @param maps A list of maps to be included. E.g. the output of `load_maps`.
arrange_data <- function(maps, duplicated_srs = 'yes', dilution_stepsize = 0) {
  # Name the maps
  for (i in seq_along(maps)) mapName(maps[[i]]) <- names(maps)[i]

  # Make sure sera names are unique
  for (i in seq_along(maps)) srNames(maps[[i]]) <- paste(names(maps)[i], srNames(maps[[i]]), sep = ':')

  # Equate similar antigens
  ag_groupings <- list(
    'D614G' = c('D614G', 'B.1')
  )
  maps <- lapply(maps, groupAgs, groupings = ag_groupings)

  if (duplicated_srs == 'yes') {
  # Include only serum groups present in at least 2 maps
    all_sr_groups <- as.character(unlist(lapply(maps, \(map) unique(srGroups(map)))))
    duplicated_sr_groups <- unique(all_sr_groups[duplicated(all_sr_groups)])
    maps <- lapply(maps, \(map) subsetMap(map, sera = srGroups(map) %in% duplicated_sr_groups))
  }

  # Get antigens that appear in at least 2 maps
  all_ags <- as.character(unlist(lapply(maps, \(map) unique(agNames(map)))))
  duplicated_ags <- unique(all_ags[duplicated(all_ags)])

  # Merge the maps
  merged_map <- mergeMaps(maps)
  # Overwrite the dilutionstepsize of 1 that mergeMaps assigns when merging maps with different dilution stepsizes.
  dilutionStepsize(merged_map) <- dilution_stepsize

  merged_map

}
