# Here I have copied a lot of source code from the gg3D website to give me more
# graphics control.

#' @importFrom magrittr "%>%"
#' @noRd
#' @references
#' \url{https://github.com/AckerDWM/gg3D/blob/master/R/labs_3D.R}
Label3D <- ggplot2::ggproto("Label3D", ggplot2::Stat,
        compute_group = function(data, scales, theta=135, phi=60, 
                            labs=c("x-axis", "y-axis", "z-axis"),
                            Xend=NA, Yend=NA, Zend=NA, common_length=1
                        ) {
                if (length(common_length)==1){Xend <- common_length
                Yend <- common_length; Zend <- common_length
                }else if (length (common_length) ==3){
                        Xend <- common_length[1]
                        Yend <- common_length[2]
                        Zend <- common_length[3]
                }
                pmat <- plot3D::perspbox(z=diag(2), plot=F, theta=theta, phi=phi)
                XY <- plot3D::trans3D(
                        x = c(Xend,0,0),
                        y = c(0,Yend,0),
                        z = c(0,0,Zend),
                pmat = pmat) %>% data.frame() %>% dplyr::mutate (label=labs)
        },
        required_aes = c("x", "y", "z")
)

#' 3D Axis Labels
#'
#' This function adds 3D axis labels to ggplot2 plots.
#'
#' @param theta The azimuthal direction in degrees.
#' @param phi The colatitude in degrees.
#' @param labs The labels to add. A vector of three where the first
#' element is x, the second is y, and the third is z.
#' @param ... Arguements passed on to layer.
#' @param common_length the proportion that the arrow axis occupies on the
#' entire axis
#' These are often aesthetics, used to set an
#' aesthetic to a fixed value, like color = "red" or size = 3.
#' @importFrom ggplot2 aes
Lab3D <- function(mapping = aes(group=1), data = NULL, geom = "text",
                    position = "identity", na.rm = FALSE, show.legend = NA,
                    inherit.aes = TRUE, AP=NULL, ...) {
        AP <- return_aes_param (AP)
        ggplot2::layer(
                stat = Label3D, data = data, mapping = mapping, geom = geom,
                position = position, show.legend = FALSE, inherit.aes = inherit.aes,
                params = list(na.rm = na.rm, size=AP$point_fontsize, ...)
        )
}

#' @importFrom magrittr "%>%"
Stat3D_proto <- ggplot2::ggproto("Stat3D", ggplot2::Stat,
        setup_params = function(data, params) {
                params$xrange <- range(data$x)
                params$yrange <- range(data$y)
                params$zrange <- range(data$z)
                params
                },
                compute_group = function(data, scales, theta=135, phi=60,
                                         xrange=c(0,1), yrange=c(0,1),
                                         zrange=c(0,1)) {
                data = data %>%
                dplyr::mutate(
                x = scales::rescale(x, from=xrange, to=c(0,1)),
                y = scales::rescale(y, from=yrange, to=c(0,1)),
                z = scales::rescale(z, from=zrange, to=c(0,1)))
                pmat = plot3D::perspbox(z=diag(2), plot=F, theta=theta, phi=phi)
                XY = plot3D::trans3D(
                        x = data$x,
                        y = data$y,
                        z = data$z,
                        pmat = pmat) %>% data.frame()
                        data$x = XY$x
                        data$y = XY$y
                        return (data)
                },
        required_aes = c("x", "y", "z")
)

#' Draw 3D Geoms
#'
#' This function adds 3D geoms such as points and paths to a ggplot2 plot.
#'
#' @param theta The azimuthal direction in degrees.
#' @param phi The colatitude in degrees.
#' @param geom The geom type to use *ie. "point", "path", "line"*
#' @param ... Arguements passed on to layer.
#' @references
#' \url{https://github.com/AckerDWM/gg3D/blob/master/R/stat_3D.R}
Stat3D <- function(mapping = NULL, data = NULL, geom = "point",
                    position = "identity", na.rm = FALSE, show.legend = NA,
                    inherit.aes = TRUE, ...) {
        ggplot2::layer(
                stat = Stat3D_proto, data = data, mapping = mapping, geom = geom,
                position = position, show.legend = show.legend, inherit.aes = inherit.aes,
                params = list(na.rm = na.rm, ...)
        )
}
