# This script cleans the immunofluorescence quantification data and combines
# them all into a single dataframe
library (magrittr)
devtools::load_all ('../..', export_all=F)
root<- '/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/results/imaging'

# where the summary data file is located
all_files <- list.files (path=paste (root, '/data/proc_data', sep=''),
                         pattern='*_sum.csv', full.names=T)
all_files %>% as.list () %>% lapply (read.csv) %>% do.call(what=rbind) -> all_data

# quality control
all_data$nuc_circ <- 4*pi*all_data$area/all_data$perimeter^2
features <- c('cell_volume', 'major_axis_length', 'minor_axis_length',
              'structure_volume', 'nuc_circ')
#quantile_plot (all_data, features, y_log=T) 
dim (all_data)
# [1] 42832    32

# append the condition labels
# I have been single blinded throughout the quantifications until this stage
cond_name <- read.csv (paste (root, 'signaling_blinding.csv', sep='/'))
series_ID <- paste (cond_name$experiment, cond_name$series_num, sep='_')
data_ID <- paste (all_data$condition, all_data$series, sep='_')
all_data$Date <- all_data$condition
all_data$condition <- cond_name$condition [match (data_ID, series_ID)]

keep_columns <- c('condition', 'Date', 'major_axis_length.struct',
                  'minor_axis_length.struct', 'perimeter.struct',
                  'area.struct', 'major_axis_length', 'minor_axis_length',
                  'perimeter', 'area', 'series', 'HLAG', 'CGB', 'TFAP2C',
                  'nuc_circ', 'HLAG_pos', 'CGB_pos', 'TFAP2C_pos', 'experimenter')
for (genes in c('HLAG', 'CGB', 'TFAP2C')){
        all_data [, paste (genes, 'pos', sep='_')] <- all_data [, genes] > 15
}
all_data %>% dplyr::mutate (experimenter= ifelse (grepl ('adam', Date), 'Adam', 'Dylan')) -> all_data
all_data %>% dplyr::select (dplyr::all_of (keep_columns)) %>% 
        write.csv (paste (paste (root, '/data/proc_data/all_cleaned2.csv', sep='')))
