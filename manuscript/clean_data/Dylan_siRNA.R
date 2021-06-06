library (magrittr)
root <- '/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/results/imaging/clono'
# combine all data together
all_files <- list.files (root, recursive=T)
all_files <- grep ('^.*/proc_data/all_sum.csv$', all_files, value=T)

all_lists <- list ()
for (i in 1:length(all_files)){
        all_lists [[i]] <- data.table::fread (paste (root, all_files[i], 
                                                     sep='/')) %>% data.frame ()
        # extract experimental condition from the file name
        all_lists [[i]]$Date <- gsub ('/proc_data/all_sum.csv', '', all_files[i])
}

all_data <- do.call (rbind, all_lists)
dim (all_data)
# [1] 2790462      14

# need to filter out the GFP trial in the first run
all_data %>% dplyr::filter (!(Date=='clono' & condition=='gfp')) -> all_data
replace_vec <- c('clono'='run1', 'gfp_RUN1_REDO'='run1', 'clono2'='run2',
                 'clono3'='run3', 'clono4'='run4', 'clono5'='run5',
                 'clono6'='run6', 'clono7'='run5', 'clono8'='run7',
                 'clono9'='run8')

all_data$Date <- replace_vec [match (all_data$Date, names (replace_vec))]
all_data %>% dplyr::mutate (condition = toupper (condition)) %>%
        dplyr::select (!V1) -> all_data

# some condition names have been mis-spelt
all_data$condition [all_data$condition=='NFE2L'] <- 'NFE2L3'
all_data$condition [all_data$condition=='PTX1'] <- 'PITX1'

# run5 and run6 share the same GFP control
gfp5 <- all_data[all_data$condition=='GFP' & all_data$Date =='run5',] 
gfp5$Date <- 'run6'
gfp7 <- all_data[all_data$condition=='GFP' & all_data$Date =='run7',] 
gfp7$Date <- 'run8'
all_data <- do.call(rbind, list (all_data, gfp5, gfp7))
write.csv (all_data, paste (root, 'combined_nuclei.csv', sep='/'))
