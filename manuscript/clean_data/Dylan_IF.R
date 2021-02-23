# This script cleans the immunofluorescence quantification data and combines
# them all into a single dataframe
library (magrittr)
devtools::load_all ('../..', export_all=F)
root<- '/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/results/imaging'
# where the summary data file is located
all_data <- read.csv (paste (root, '/data/proc_data/all_sum.csv', sep=''))

# quality control
all_data$nuc_circ <- 4*pi*all_data$area/all_data$perimeter^2
features <- c('cell_volume', 'major_axis_length', 'minor_axis_length',
              'structure_volume', 'nuc_circ')
quantile_plot (all_data, features, y_log=T) 
dim (all_data)
# [1] 31934    32

# check the QC after filtering
all_data2 <- all_data [all_data$nuc_circ > 0.3,]
quantile_plot (all_data2, features, y_log=T) 
dim (all_data2)
# [1] 31836    32

# append the condition labels
# I have been single blinded throughout the quantifications until this stage
cond_name <- read.csv (paste (root, 'signaling_blinding.csv', sep='/'))
all_data2$condition <- cond_name$condition [match (all_data2$series, cond_name$series_num)]
keep_columns <- c('condition', 'major_axis_length.struct',
                  'minor_axis_length.struct', 'perimeter.struct',
                  'area.struct', 'major_axis_length', 'minor_axis_length',
                  'perimeter', 'area', 'series', 'HLAG', 'CGB', 'TFAP2C',
                  'nuc_circ', 'HLAG_pos', 'CGB_pos', 'TFAP2C_pos')
for (genes in c('HLAG', 'CGB', 'TFAP2C')){
        all_data2 [, paste (genes, 'pos', sep='_')] <- all_data2 [, genes] > 7
}
all_data2 %>% dplyr::select (dplyr::all_of (keep_columns)) %>% 
        write.csv (paste (paste (root, '/data/proc_data/all_cleaned.csv', sep='')))
