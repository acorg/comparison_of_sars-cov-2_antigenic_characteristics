#' Code for plotting a map with highlighted and de-emphasized antigens
#' @param input_map The acmap to plot
#' @param grey The antigens to de-emphasize
#' @param highlight The antigens to highlight
#' @param alphaF The alpha of the de-emphasized antigens
#' @param sr_outline_width Width of the serum outline
#' @param xlim X-limits
#' @param ylim Y-limits
#' @param cex Point scale
#' @param prucrustes Whether procurstes arrows should be shown

greyed_map <- function(input_map, grey, highlight, alphaF,
  sr_outline_width = 0.2, xlim = c(-1, 9), ylim = c(-2, 3), cex=1.5,
  procrustes = NULL) {
  # Plot a map where important antigens are highlighted

  map <- input_map

  srOutlineWidth(map) <- sr_outline_width

  unwanted_ags <- agNames(map) %in% grey
  highlighted_ags <- agNames(map) %in% highlight
  oas <- agNames(map)[! agNames(map) %in% c(grey, highlight)]
  other_ags <- agNames(map) %in% oas

  # Style the antigens that shouldn't be highlighted or greyed out
  agFill(map)[other_ags] <- adjustcolor(agFill(map)[other_ags], alpha.f = 0.9)
  agOutline(map)[other_ags] <- adjustcolor(agOutline(map)[other_ags], alpha.f = 0.9)

  # Lower the visibility of unwanted antigens
  agFill(map)[unwanted_ags] <- adjustcolor(agFill(map)[unwanted_ags], alpha.f = alphaF)
  agOutline(map)[unwanted_ags] <- adjustcolor(agOutline(map)[unwanted_ags], alpha.f = alphaF)

  # Increase the visibility of the antigens to be highlighted
  agOutlineWidth(map)[highlighted_ags] <- 3
  agFill(map)[highlighted_ags] <- adjustcolor(agFill(map)[highlighted_ags], alpha.f = 1)

  # If procrustes arrows should be drawn, only draw them for the non-unwanted antigens.
  if (!(is.null(procrustes))) {
    procrustes_colors <- rep('black', numPoints(map))
    procrustes_colors[c(unwanted_ags)] <- adjustcolor('black', alpha.f = 0.2)
  } else {
    procrustes_colors <- rep('black', numPoints(map))
  }

  # Plot the map
  par(mar = c(0.3, 0, 0.3, 0))
  plot(map, xlim = xlim, ylim = ylim, fill.alpha = 1, plot_labels = FALSE,
       outline.alpha = 1, grid.col = '#cfcfcf', plot_stress = FALSE,
       grid.margin.col = 'black', show_error_lines = FALSE, cex=cex,
       procrustes.col = procrustes_colors, optimization_number = 1)
}


