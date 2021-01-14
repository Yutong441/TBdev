# ----------colors----------
# create matching color
setwd('..')
devtools::load_all ()
data (color_scheme)
data (orders)
color_vec <- paste ('#', as.character (color_scheme$HEX), sep='')
names (color_vec) <- color_scheme$celltype
color_vec <- color_vec [ !(names (color_vec) == '') ]

order_color <- partial_relevel (names (color_vec), c(orders$cell_order, orders$time_order))
order_color <- partial_relevel (order_color, orders$branch_order)
color_vec <- color_vec [order (order_color)]
usethis::use_data (color_vec, overwrite=T)
