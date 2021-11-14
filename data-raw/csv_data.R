# Process all csv data
# the working directory is already data-raw
data_dir <- 'config'
cell_cycle <- read.csv (paste (data_dir, 'cell_cycle.csv', sep='/'))
vivo_clust <- read.csv (paste (data_dir, 'cluster_invivo.csv', sep='/'))
color_scheme <- read.csv (paste (data_dir, 'expanded_color.csv', sep='/'))
date2CS <- read.csv (paste (data_dir, 'date2CS.csv', sep='/'))
KeggID <- read.csv (paste (data_dir, 'KEGG_ID.csv', sep='/'))
convert_name <- read.csv (paste (data_dir, 'name_conversion.csv', sep='/'))
GOsimp <- read.csv (paste (data_dir, 'simplify_GO_terms.csv', sep='/'))
usethis::use_data (cell_cycle, vivo_clust, color_scheme, date2CS, KeggID,
                   convert_name, GOsimp, overwrite=T)
