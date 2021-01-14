root <- '/mnt/c/users/Yutong/Documents/bioinformatics/reproduction/'
#devtools::create (paste (root,'TBdev', sep=''))
setwd ('..')

data (lineage_markers)
roxygen2::roxygenise()
devtools::load_all ()

root_dir <- '/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/'
root2 <- paste (root_dir, 'results/', sep='/')
merge_dir <- paste (root2, 'XLYBPZ_Dylan_dir', sep='/')
x <- load (paste (merge_dir, 'final_merged_tb.Robj', sep='/') )
all_data <- get (x)

sup_save_dir <- paste (root2, 'manuscript/figureS2', sep='/')
BRGP_prob <- read.csv (paste (sup_save_dir, 'result/BRGP_likelihood.csv', sep='/') )
show_meta <- all_data@meta.data [!is.na (all_data$MGP_PT) & !all_data$broad_type %in% c('EPI', 'PE'),]
select_cells <- c('eICM', 'aICM', 'eTB', 'iTB', 'aTB', 'hTSC_OKAE', 'hTSC_TURCO', 'cleavage',
                  'eEVT', 'aEVT', 'eCTB', 'vCTB1', 'vCTB2', 'vCTB3',
                  'hESC', 'hESC_YAN')

devtools::load_all ()
data (CT)
plot_prob_line (BRGP_prob, select_cells, CT$in_vitro_cells,
                meta=show_meta[show_meta$broad_type!='STB',], 
                vjust=1e-4, thickness=1e-4, normalize_data=T)

select_cells <- c('cleavage', 'eICM', 'aICM', 'eTB', 'iTB', 'aTB', 'hTSC_OKAE', 'hTSC_TURCO',
                  'hESC', 'hESC_YAN', 'eCTB', 'sCTB', 'aSTB1', 'aSTB2', 'aSTB3')
plot_prob_line (BRGP_prob, select_cells, CT$in_vitro_cells,
                meta=show_meta [show_meta$broad_type!='EVT',], 
                vjust=1e-4, thickness=1e-4, sel_branch='branch2', normalize_data=T)
