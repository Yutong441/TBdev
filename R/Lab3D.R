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
#' @param label The labels to add. A vector of three where the first
#' element is x, the second is y, and the third is z.
#' @param ... Arguements passed on to layer.
#' @param common_length the proportion that the arrow axis occupies on the
#' entire axis
#' These are often aesthetics, used to set an
#' aesthetic to a fixed value, like color = "red" or size = 3.
#' @importFrom ggplot2 aes
#' @export
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
