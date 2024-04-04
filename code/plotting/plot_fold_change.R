
#' Convert log2 to fold change
foldchange <- \(x) {

  xabs <- abs(x)
  foldchange <- 2^xabs
  foldchange[x < 0] <- -foldchange[x < 0]
  as.character(round(foldchange, 1))

}


#' Code to make the fold change bar plots shown in figures S36-38
#' @param data The foldchange data to plot
#' @param cols The colors
#' @param title The plot title
#' @param shapes The shapes for each dataset
#' @param highlighted_map The names of any maps to highlight by using a bigger
#' symbol
make_single_foldchange_plot_no_gmt <- function(data, cols = map_colors,
  title_ = '', shapes = map_shapes, highlighted_map = c()) {

  size_values <- rep(0.7, length(unique(data$map)))
  names(size_values) <- unique(data$map)

  for (m in highlighted_map) {
    size_values[m] <- 0.3
  }

  data %>%
    group_by(
      sr_group, ag_name
    ) %>%
    summarise(
      mean_diff = mean(mean_diff, na.rm = T)
    ) -> mean_fold_change

  data %>%
    filter(map != 'estimated gmt') %>%
    ggplot(
      aes(
        x = ag_name,
        y = mean_diff,
        color = map
      )
    ) +
    geom_hline(
      yintercept = 0,
      lty = 'dashed',
      color = 'grey'
    ) +
    geom_barcol(
      data = mean_fold_change,
      aes(
        x = ag_name,
        y = mean_diff
      ),
      fill = 'grey90',
      color = 'grey90',
      alpha = 1,
      base = 0
    ) +
    geom_pointrange(
      aes(
        ymin = mean_diff_lower,
        ymax = mean_diff_upper,
        shape = map,
        size = map
      ),
      position = position_dodge2(
        width = 0.8
      )
    ) +
    coord_cartesian(
      ylim = c(-12, 3)
    ) +
    scale_y_continuous(
      breaks = function(x){
        ceiling(x[1]):floor(x[2])
      },
      labels = rev(c('', '4', '', '1', '', '-4', '', '-16', '', '-64', '', '-256', '', '-1024', '', '-4096'))
    ) +
    theme_bw() +
    geom_vline(
      xintercept = seq(-0.5, 20.5, by = 1),
      color = 'gray80'
    ) +
    scale_x_discrete(
      limits =  rev(mean_fold_change$ag_name[order(mean_fold_change$mean_diff, na.last = NA)])
    ) +
    theme( # remove the vertical grid lines
      panel.grid.major.x = element_blank()
    ) +
    theme_mapplot() +
    scale_color_manual(
      name = '',
      values = cols
    ) +
    scale_size_manual(
      name = '',
      values = size_values
    ) +
    scale_shape_manual(
      values = shapes,
      name = '',
      labels = unlist(lapply(animal_dataset_order, function(i) pretty_labels[i]),
                                                           use.names=TRUE)
    ) +
    guides(size = 'none') +
    labs(
      y = 'Fold change',
      x = '',
      title = title_
    ) +
    theme(
      axis.text.x = element_text(size = 20),
      axis.title = element_text(size = 20),
      axis.text.y = element_text(size=20),
      legend.text = element_text(size=20),
      plot.title = element_text(size = 25)
    ) -> gp
}


#' Function to plot fold changes grouped by map and with added slopes.
#' Shown in figures S34-35 of the SOM.
#' @param data The fold change data to plot
#' @param slope_data The slope data
#' @param labels X axis labels
#' @param line_indexes The indices of the vertical lines separating each map
fold_change_by_map_with_slopes <- function(
    data, slope_data, labels, line_indexes) {
  data %>%
  ggplot(
  ) +
  geom_line(
    data = slope_data,
    aes(
      x = x_loc,
      y = foldchange
    ),
    color = 'black',
    alpha = 0.9,
    size = 1
  ) +
  geom_hline(
    yintercept = 0,
    color = 'grey',
    linetype = 'dashed'
  ) +
  geom_pointrange(
    data = data,
    aes(
      x = x_loc,
      y = mean_diff,
      ymin = mean_diff_lower,
      ymax = mean_diff_upper,
      color = ag_name,
      shape = ag_name
    ),
    #shape = 16,
    size = 1.4
  ) +
  scale_color_manual(
    values = ag_colors,
    name = 'Variant'
  ) +
  scale_shape_manual(
    values = ag_shapes,
    name = 'Variant'
  ) +
  theme_mapplot() +
  geom_vline(
    xintercept = c(seq(0.975, 16.975, by = 1),
    color = 'gray30'
  ) +
  geom_hline(
    yintercept = 0,
    color = 'black'
  ) +
  labs(
    x = '',
    y = 'Fold change'
  ) +
  theme_bw() +
  theme( # remove the vertical grid lines
    panel.grid.major.x = element_blank(),
    panel.border = ggplot2::element_rect(color = 'black', fill = NA, size = 1)
  ) +
  theme(text = element_text(size = 35)) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5)) +
  theme(panel.spacing.y = unit(2, 'lines')) +
  scale_x_continuous(
    breaks = seq(0.475, 17.475, by = 1),
    labels = labels
  ) +
  geom_vline(
    xintercept = line_indexes,
    color = 'gray30',
    size = 3
  ) +
  coord_cartesian(
    ylim = c(-9, 2.5),
    xlim = c(0.79, 17.16)
  ) +
  scale_y_continuous(
    breaks = function(x) ceiling(min(x)):floor(max(x)),
    labels = log2foldchange
  ) +
  theme(
    axis.text.x = element_markdown(vjust = 0.9, hjust = 0.5, size = 30),
    strip.background = element_rect(colour='white', fill='white'),
    strip.text = element_text(size = 37),
    legend.key.size = unit(1.2, 'cm')
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 2), ncol = 1),
    shape = guide_legend(override.aes = list(size = 2), ncol = 1)
  ) +
  facet_rep_wrap(~sr_group, ncol = 1, repeat.tick.labels = T) -> gp
}


#' Function to plot fold changes grouped by map and with added slopes.
#' Shown in figure 2.
#' @param data The fold change data to plot
#' @param slope_data The slope data
#' @param slope_data_ind The slope data
#' @param labels X axis labels
#' @param line_indexes The indices of the vertical lines separating each map
fold_change_by_map_with_slopes_main <- function(
    data, slope_data, slope_data_ind, labels, line_indexes) {
  data %>%
  ggplot(
  ) +
  geom_line(
    data = slope_data,
    aes(
      x = x_loc,
      y = foldchange
    ),
    color = 'black',
    size = 1.5
  ) +
  geom_line(
  data = slope_data_ind,
  aes(
    x = x_loc,
    y = foldchange
  ),
  color = 'gray60',
  size = 1
  ) +
  geom_hline(
    yintercept = 0,
    color = 'grey',
    linetype = 'dashed'
  ) +
  geom_pointrange(
    data = data,
    aes(
      x = x_loc,
      y = mean_diff,
      ymin = mean_diff_lower,
      ymax = mean_diff_upper,
      color = ag_name,
      shape = ag_name
    ),
    size = 1.4
  ) +
  scale_color_manual(
    values = ag_colors,
    name = 'Variant'
  ) +
  scale_shape_manual(
    values = ag_shapes,
    name = 'Variant'
  ) +
  theme_mapplot() +
  geom_vline(
    xintercept = c(seq(0.975, 16.975, by = 1),
    color = 'gray30'
  ) +
  geom_hline(
    yintercept = 0,
    color = 'black'
  ) +
  labs(
    x = '',
    y = 'Fold change'
  ) +
  theme_bw() +
  theme( # remove the vertical grid lines
    panel.grid.major.x = element_blank(),
    panel.border = ggplot2::element_rect(color = 'black', fill = NA, size = 1)
  ) +
  theme(text = element_text(size = 35)) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5)) +
  theme(panel.spacing.y = unit(2, 'lines')) +
  scale_x_continuous(
    breaks = seq(0.475, 17.475, by = 1),
    labels = labels
  ) +
  geom_vline(
    xintercept = line_indexes,
    color = 'gray30',
    size = 3
  ) +
  coord_cartesian(
    ylim = c(-9, 2.5),
    xlim = c(0.79, 17.16)
  ) +
  scale_y_continuous(
    breaks = function(x) ceiling(min(x)):floor(max(x)),
    labels = log2foldchange
  ) +
  theme(
    axis.text.x = element_markdown(vjust = 0, hjust = 0.5, size = 35),
    strip.background = element_rect(colour='white', fill='white'),
    strip.text = element_text(size = 37),
    legend.key.size = unit(1.2, 'cm')
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 2), ncol = 1),
    shape = guide_legend(override.aes = list(size = 2), ncol = 1)
  ) +
  facet_wrap(~sr_group, ncol = 1) -> gp
}

