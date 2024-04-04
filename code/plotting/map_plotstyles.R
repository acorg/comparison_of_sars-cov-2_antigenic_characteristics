
#' @export
theme_mapplot <- function(
  axis.text.x.angle = 45,
  axis.text.x.hjust = 1,
  axis.text.x.vjust = 1
  ) {
  theme(
    panel.background = ggplot2::element_blank(),
    plot.background = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(color = 'grey50', fill = NA),
    axis.text.x = ggplot2::element_text(
      angle = axis.text.x.angle,
      hjust = axis.text.x.hjust,
      vjust = axis.text.x.vjust
    ),
    panel.grid = element_line(
      colour = '#eeeeee'
    ),
    axis.line = ggplot2::element_line(
      colour = 'grey80'
    )
  )
}


#' @export
scale_y_titer <- function(threshold = '<10', axisname = 'Titer') {

  logthreshold <- Racmacs:::log_titers(threshold, 1)

  scale_y_continuous(
    name = axisname,
    breaks = function(x) {
      floor(logthreshold):ceiling(max(x))
    },
    labels = function(x) {
      output <- 2^x*10
      output[x == floor(logthreshold)] <- threshold
      output
    },
    minor_breaks = NULL
  )

}
