# basic settings
# !!! must be run after running the 'color_vec.R' file and 'cell_markers.R'
#highlight_font <- gpar (fontface='bold', fontsize=18)
setwd ('..')
data (color_vec)
data (orders)

remove_keys <- c('disease', 'infection', 'shigellosis', 'atherosclerosis',
                 'measles', 'bacteria', 'pertussis', 'sclerosis', 'hepatitis',
                 'ataxia', 'syndrome', 'addiction', 'Toxoplasmosis',
                 'Leishmaniasis', 'diabetes', 'tuberculosis', 'myocarditis',
                 'Amoebiasis', 'Asthma', 'influenza', 'legionellosis', 'lupus',
                 'malaria', 'arthritis', 'HIV')

format_conf <- list (
        fontsize = 15,
        point_fontsize = 18/3, #for geom_text
        pointsize = 3,
        legend_point_size = 5,
        normal_shape = 21, # circle, has to be larger than 20
        point_edge_color = 'gray',
        edge_stroke= 0.8,
        highlight_shape = 24, # triangle
        font_fam = 'Arial',
        # for grid::gpar argument, 'sans' mans 'Arial'
        gfont_fam = 'sans', 
        highlight_font = list(fontface='bold', fontsize=18),
        ridge_alpha=0.5,

        # color settings
        heatmap_color = c("#00B0F0", "#FFF2F1", "#FF0000"), #blue, white, red
        color_vec = color_vec,
        date_color_vec = c(in_vitro='#FF0000', no_dates='#FF0000'),
        palette = 'viridis', # for `ggplot2::scale_color_continuous`

        remove_keys = remove_keys,
        # arrow settings
        arrow_angle=15,
        arrow_length=0.05, 
        arrow_length_unit='npc',
        arrow_type = 'closed',
        arrow_thickness= 0.5,
        arrow_font_amp=0.8,
        arrow_linejoin = 'mitre' #sharp border
)
format_conf <- append (format_conf, orders)
usethis::use_data (format_conf, overwrite=T)
