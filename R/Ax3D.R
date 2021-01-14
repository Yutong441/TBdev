#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @importFrom plot3D trans3D
#' @references
#' \url{https://github.com/AckerDWM/gg3D/blob/master/R/axes_3D.R}
#' @author AckerDWM (modified by Yutong Chen)
coord_transform <- function (theta, phi, Xend, Yend, Zend){
    pmat <- plot3D::perspbox(z=diag(2), plot=F, theta=theta, phi=phi)
    x_axis = trans3D(x = c(0, Xend), y = 0, z = 0, pmat = pmat) %>%
      data.frame() %>% mutate(axis="x") %>% mutate (pos=c('start', 'end')) 
    y_axis = trans3D(x = 0, y = c(0, Yend), z = 0, pmat = pmat) %>%
      data.frame() %>% mutate(axis="y")%>% mutate (pos=c('start', 'end'))
    z_axis = trans3D(x = 0, y = 0, z = c(0, Zend), pmat = pmat) %>%
      data.frame() %>% mutate(axis="z")%>% mutate (pos=c('start', 'end'))
    return (dplyr::bind_rows(x_axis, y_axis, z_axis))
}

Axes3D <- ggplot2::ggproto("Axes3D", ggplot2::Stat,
        compute_group = function(data, scales, theta=135, phi=60,
                            Xend=NA, Yend=NA, Zend=NA, common_length=1) {
                if (!is.na(common_length)){Xend <- common_length
                Yend <- common_length; Zend <- common_length
                }
                return (coord_transform (theta, phi, Xend, Yend, Zend) )
        },
        required_aes = c("x", "y", "z")
)

#' Draw 3D Axes
#'
#' This function adds 3D axes to a ggplot2 plot.
#'
#' @param theta The azimuthal direction in degrees.
#' @param phi The colatitude in degrees.
#' @param ... Arguements passed on to layer.
#' These are often aesthetics, used to set an
#' aesthetic to a fixed value, like color = "red" or size = 3.
#' @importFrom ggplot2 aes
#' @export
Ax3D = function(mapping = aes(group=1), data = NULL, geom = "path",
                    position = "identity", na.rm = FALSE, show.legend = NA,
                    inherit.aes = TRUE, ...) {
        ggplot2::layer(
    stat = Axes3D, data = data, mapping = mapping, geom = geom,
    position = position, show.legend = FALSE, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}

#' @importFrom magrittr %>%
#' @noRd
Segment3D <- ggplot2::ggproto("Axes3D", ggplot2::Stat,
        compute_group = function(data, scales, theta=135, phi=60,
                            Xend=NA, Yend=NA, Zend=NA,
                            common_length=1) {
                if (!is.na(common_length)){Xend <- common_length
                Yend <- common_length; Zend <- common_length
                }
                trans <- coord_transform (theta, phi, Xend, Yend, Zend)
                trans %>% dplyr::filter (pos == 'start') -> start_trans
                trans %>% dplyr::filter (pos == 'end') %>% dplyr::select (!axis) %>%
                        magrittr::set_colnames(c('xend', 'yend', 'axis'))-> end_trans
                com_trans <- cbind (start_trans, end_trans)
                return (com_trans)
        },
        required_aes = c("x", "y", "z")
)

#' Extention to `Ax3D` with arrow heads
#'
#' @importFrom ggplot2 aes
#' @export
Seg3D <- function(mapping = aes(group=1), data = NULL, geom = "segment",
                  position = "identity", na.rm = FALSE, show.legend = NA,
                  inherit.aes = TRUE, AP=NULL, ...) {
        AP <- return_aes_param (AP)
        ggplot2::layer(
                stat = Segment3D, data = data, mapping = mapping, geom = geom,
                position = position, show.legend = FALSE, inherit.aes = inherit.aes,
                params = list(na.rm = na.rm, arrow=get_arrow(AP), size=AP$arrow_thickness, 
                              linejoin=AP$arrow_linejoin, ...)
        )
}
