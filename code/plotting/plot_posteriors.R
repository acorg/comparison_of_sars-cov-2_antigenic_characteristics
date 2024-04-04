#' Plot a violin plot of the samples from the posterior distributions
#' @param ci_data: The calculated confidence intervals, for plotting the confidence intervals
#' @param draws_data: The raw samples for plotting the violin plots.
#' @param map_order: The order in which the datasets should be plotted.
#' @param map_colors: The colors in which the datasets should be plotted.
#' @param xlabel: The label for the x-axis label.
plot_posteriors_hpdi_violin <- function(ci_data, draws_data, map_order,
  map_cols, xlabel, xlim = c(-10, 10)
  ) {
  ggplot() +
    geom_violin(
      data = draws_data,
      aes(
        y = map,
        x = value,
        fill = map
      ),
      color = 'black',
      alpha = 0.5,
      size = 0.3,
      scale = 'count'
    ) +
    geom_linerange(
      data = ci_data,
      aes(
        y = map,
        xmin = ci_low,
        xmax = ci_high,
        color = map
      ),
      size = 1.5,
      alpha = 1
    ) +
    geom_point(
      data = ci_data,
      aes(
        y = map,
        x = mean_,
        color = map
      ),
      size = 2,
      alpha = 1
    ) +
    geom_point(
      data = ci_data,
      aes(
        y = map,
        x = mean_,
      ),
      size = 1,
      alpha = 0.5,
      fill = 'white',
      color = 'white'
    ) +
    geom_vline(
      xintercept = 0,
      linetype = 'dashed',
      color = 'black',
      alpha = 0.5
    ) +
    scale_y_discrete(
      limits = map_order
    ) +
    labs(
      x = xlabel,
      y = ''
    ) +
    scale_color_manual(
      values = map_cols
    ) +
    scale_fill_manual(
      values = map_cols
    ) +
    coord_cartesian(
      xlim = xlim
    ) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 15),
      axis.text.x = element_text(size = 12),
      axis.title.x = element_text(size = 15)
    ) +
    guides(color = 'none', size = 'none', fill = 'none') -> gp
}

