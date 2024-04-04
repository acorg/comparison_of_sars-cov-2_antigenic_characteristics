
#' @export
GeomBarCol <- ggplot2::ggproto(
    'GeomBarCol', ggplot2::GeomCol,
    draw_panel = function(self, data, panel_params, coord, width = NULL, base = NULL, flipped_aes = FALSE) {

        # Make edits to data
        y <- data$ymax
        y[data$ymin < 0] <- data$ymin[data$ymin < 0]
        data$ymin <- pmin(y, base)
        data$ymax <- pmax(y, base)

        # Pass to the geom to plot
        ggplot2::ggproto_parent(ggplot2::GeomRect, self)$draw_panel(data, panel_params, coord)

    }
)

#' @export
geom_barcol <- function (
    mapping = NULL,
    data = NULL,
    position = 'stack',
    ...,
    width = NULL,
    base = 0,
    na.rm = FALSE,
    show.legend = NA,
    inherit.aes = TRUE
) {

    ggplot2::layer(
        data = data,
        mapping = mapping,
        stat = 'identity',
        geom = GeomBarCol,
        position = position,
        show.legend = show.legend,
        inherit.aes = inherit.aes,
        params = list(
            width = width,
            na.rm = na.rm,
            base = base,
            ...
        )
    )

}

