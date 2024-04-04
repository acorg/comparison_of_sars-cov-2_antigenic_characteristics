#' Code to plot titers used in Figure 1
#' @param sr_group_plotdata GMTs per map, serum group and antigen
#' @param bar_data GMTs per serum group and antigen
#' @param colors Colors for each dataset
#' @param shapes Shape for each dataset
#' @param included_maps Maps to be included in the plot
#' @param highlighted_map Map(s) to highlight by a bigger symbol
#' @param ag_order The order of the antigens
#' @param bar_col The color of the bars
#' @param bar_alpha The opacity of the bars
plot_adjusted_titer_comparison_gmt <- function(
  sr_group_plotdata,
  bar_data,
  colors = map_colors,
  shapes = map_shapes,
  included_maps = unique(sr_group_plotdata$map),
  highlighted_map = c(),
  ag_order = unique(sr_group_plotdata$ag_name),
  bar_color = 'grey90',
  bar_alpha = 1
  ) {

  ag_order <- ag_order[(ag_order %in% unique(sr_group_plotdata$ag_name))]

  size_values <- rep(0.7, length(unique(sr_group_plotdata$map)))
  names(size_values) <- unique(sr_group_plotdata$map)

  for (m in highlighted_map) {
    size_values[m] <- 0.3
  }

  bar_data %>%
    filter(
      paste(sr_group, ag_name) %in% paste(sr_group_plotdata$sr_group, sr_group_plotdata$ag_name)
    ) %>%
    subset(estimate != 7) -> bar_data

  bar_data$ag_name <- factor(bar_data$ag_name, levels = ag_order)
  sr_group_plotdata$ag_name <- factor(sr_group_plotdata$ag_name, levels = ag_order)

  sr_group_plotdata %>%
    left_join(bar_data, by = c('sr_group', 'ag_name')) %>%
    subset(map %in% included_maps) %>%
    filter(
      paste(sr_group, ag_name) %in% paste(bar_data$sr_group, bar_data$ag_name)
      ) %>%
    ggplot(
      aes(
        x = ag_name,
        y = gmt
      )
    ) +
    geom_barcol(
      data = bar_data,
      aes(
        x = ag_name,
        y = estimate
      ),
      fill = bar_color,
      alpha = bar_alpha,
      base = -10
    ) +
    geom_pointrange(
      aes(
        ymin = gmt_lower,
        ymax = gmt_upper,
        color = map,
        shape = map,
        size = map
      ),
      position = position_dodge2(
        width = 0.8
      )
    ) +
    scale_color_manual(
      values = colors
    ) +
    scale_size_manual(
      values = size_values
    ) +
    scale_shape_manual(
      values = shapes
    ) +
    guides(
      size = 'none',
      shape = 'none',
      fill = guide_legend(override.aes = list(alpha = 1))
    ) +
    coord_cartesian(
      ylim = c(0, 14)
    ) +
    scale_y_titer(threshold = '<5') +
    theme(
      text = element_text(size = 20),
      legend.position = 'top',
      legend.key = element_rect(fill = 'white'),
      strip.background = element_blank()
    ) +
    facet_wrap(
      ~sr_group,
      ncol = 1
    ) +
    labs(
      x = 'Antigen',
      color = ''
    ) +
    theme(
      panel.background = element_rect(fill = 'white', color = 'white')
    ) +
    theme_mapplot()

}
