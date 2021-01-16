setwd ('..')
devtools::load_all ()
library (tidyverse)
# ----------plotting settings----------
format_conf <- list (
        fontsize = 15,
        point_fontsize = 18/3, #for geom_text
        pointsize = 3,
        legend_point_size = 5,
        normal_shape = 21, # circle, has to be larger than 20
        point_edge_color = 'white',
        highlight_shape = 24, # triangle
        font_fam = 'Arial',
        # for grid::gpar argument, 'sans' mans 'Arial'
        gfont_fam = 'sans', 
        highlight_font = list(fontface='bold', fontsize=18),
        ridge_alpha=0.5,

        # color settings
        heatmap_color = c("#00B0F0", "#FFF2F1", "#FF0000"), #blue, white, red
        date_color_vec = c(in_vitro='#FF0000'),
        palette = 'viridis', # for `ggplot2::scale_color_continuous`

        # arrow settings
        arrow_angle=30,
        arrow_length=0.01, 
        arrow_length_unit='npc',
        arrow_type = 'closed',
        arrow_thickness= 1.5,
        arrow_linejoin = 'mitre' #sharp border
)

# ----------celltype settings----------
# define a list of different cell types for filtering
CT <- list(
        TB_lineage = c('STB', 'EVT', 'ICM', 'CTB', 'TE'),
        non_emb_lineage = c('STR', 'STROMA', 'RESTROMA', 'OVIDUCT', 'NEONATE',
                             'GLAND', 'REGLAND', 'F', 'NL', 'NOEXPRESSION', 'MIX', 'MYO'),
        pre_imp_lineage = c('Oocyte', 'Zy', '2C', '4C', '8C', 'cMor', 'PE',
                            'HYP', 'VE', 'EPI', 'PSA-EPI', 'preICM',
                            'cleavage'),
        in_vitro_cells = c('hESC', 'hES', 'ESC', 'hESC_YAN', 'hESC', 'maESC',
                           'hTSC', 'hTSC_OKAE', 'hTSC_TURCO')
)
CT$non_TB_lineage <- c(CT$non_emb_lineage,
                              CT$pre_imp_lineage,
                              CT$in_vitro_cells)

# ----------ordering settings----------
# regulate the appearance of cells in legends or other features
# can be used in conjunction with `partial_relevel`
all_orders <- list (
        cell_order = c('Oocyte', 'Zy', '2C', '4C', '8C', 'cMor', 'Blast',
                        'cleavage', 'ICM', 'EPI', 'PSA-EPI','PE', 'HYP', 'TB',
                        'TE', 'CTB', 'STB', 'EVT', 'STR')
)

# ----------color settings----------

# To add more colors, add 3 fields separated by commas
# The first indicates the color names (not essential)
# The second is the HEX code without '#'
# The third is the names of the features to be colored
# Using more than one HEX code for different features is allowed
color_vec <- c(
        "Light-green1",  "56E600",  "cleavage"   , 
        "Light-green1",  "56E600",  "Oocyte"     , 
        "Light-green2",  "48BF00",  "Zy"         , 
        "Light-green3",  "399900",  "2C"         , 
        "Green1"      ,  "00E639",  "4C"         , 
        "Green2"      ,  "00BF30",  "8C"         , 
        "Green3"      ,  "009926",  "cMor"       , 
        "Green3"      ,  "009926",  "Blast"      , 
        "Green3"      ,  "009926",  "eICM"       , 
        "Turquoise1"  ,  "00E6E6",  "ICM"        , 
        "Turquoise1"  ,  "00E6E6",  "aICM"       , 
        "Turquoise2"  ,  "00BFBF",  "EPI"        , 
        "Turquoise2"  ,  "00BFBF",  "EPI1"       , 
        "Turquoise2"  ,  "00BFBF",  "EPI2"       , 
        "Turquoise3"  ,  "009999",  "PSA-EPI"    , 
        "Light-blue1" ,  "0039E6",  "sCTB"       , 
        "Light-blue1" ,  "0039E6",  "sCTB1"      , 
        "Light-blue2" ,  "0233BF",  "sCTB2"      , 
        "Light-blue3" ,  "002699",  "sCTB3"      , 
        "Light-blue3" ,  "002699",  "eSTB"       , 
        "Dark-blue1"  ,  "2600E6",  "STR"        , 
        "Dark-blue1"  ,  "2600E6",  "aSTB1"      , 
        "Dark-blue2"  ,  "2000BF",  "aSTB2"      , 
        "Dark-blue3"  ,  "2C1599",  "aSTB3"      , 
        "Dark-blue4"  ,  "1A0873",  "STB"        , 
        "Dark-blue4"  ,  "1A0873",  "aSTB4"      , 
        "Dark-blue5"  ,  "13084D",  ""           , 
        "Purple1"     ,  "921FE6",  "CTB"        , 
        "Purple1"     ,  "921FE6",  "eCTB"       , 
        "Purple2"     ,  "7813BF",  ""           , 
        "Purple3"     ,  "5B0399",  ""           , 
        "Pink1"       ,  "E605A4",  "TB"         , 
        "Pink1"       ,  "E605A4",  "eTB"        , 
        "Pink1"       ,  "E605A4",  "TE"         , 
        "Pink2"       ,  "BF0489",  "iTB"        , 
        "Pink3"       ,  "99036D",  "aTB"        , 
        "Light-red1"  ,  "F04C04",  "hTSC_TURCO" , 
        "Light-red2"  ,  "BF3C04",  "hTSC_OKAE"  , 
        "Light-red2"  ,  "BF3C04",  "hTSC"       , 
        "Light-red3"  ,  "992F03",  "vCTB1"      , 
        "Red1"        ,  "E60505",  "vCTB2"      , 
        "Red1"        ,  "E60505",  "vCTB3"      , 
        "Red2"        ,  "BF0504",  "eEVT"       , 
        "Red2"        ,  "BF0504",  "EVT"        , 
        "Red3"        ,  "960303",  "aEVT"       , 
        "Orange1"     ,  "E68600",  ""           , 
        "Orange2"     ,  "BF7104",  ""           , 
        "Orange3"     ,  "995B03",  ""           , 
        "Gold1"       ,  "E6B500",  "PE"         , 
        "Gold2"       ,  "BF9600",  ""           , 
        "Gold3"       ,  "967700",  ""           , 
        "Yellow1"     ,  "E6E600",  "hESC"       , 
        "Yellow1"     ,  "E6E600",  "hESC1"      , 
        "Yellow2"     ,  "BFBF04",  "hESC_YAN"   , 
        "Yellow2"     ,  "BFBF04",  "hESC2"      , 
        "Yellow2"     ,  "BFBF04",  "hES"        , 
        "Yellow3"     ,  "999903",  ""           , 
        "Gray"        ,  "808080",  "MIX"        , 
        "Gray"        ,  "808080",  "unknown"    , 
        "Gray"        ,  "808080",  "uCTB"        
)

color_vec <- matrix (color_vec, ncol=3, byrow=T)
color_vec [,3:2] %>% data.frame () %>% deframe () -> color_vec
names (color_vec) <- partial_relevel (names (color_vec), all_orders$cell_order)

# ----------GO/KEGG/Reactome settings----------
# prevent showing certain terms in GO/KEGG/Reactome analysis
remove_keys <- c('disease', 'infection', 'shigellosis', 'atherosclerosis',
                 'measles', 'bacteria', 'pertussis', 'sclerosis', 'hepatitis',
                 'ataxia', 'syndrome', 'addiction', 'Toxoplasmosis',
                 'Leishmaniasis', 'diabetes', 'tuberculosis', 'myocarditis',
                 'Amoebiasis', 'Asthma', 'influenza', 'legionellosis', 'lupus',
                 'malaria', 'arthritis')

format_cf <- c(format_conf, CT, all_orders, list (remove_keys=remove_keys, color_vec=color_vec))
rm (all_orders, format_conf, remove_keys, color_vec, CT)
