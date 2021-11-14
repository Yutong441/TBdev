library (tidyverse)
devtools::load_all ('..', export_all=T)
root_dir <- '/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/results'
root <- paste (root_dir, 'imaging/clono', sep='/')
save_dir <- paste (root_dir, 'manuscript/figure4/clonogenicity', sep='/')
sup_save_dir <- paste (root_dir, 'manuscript/figureS5', sep='/')
all_data <- data.table::fread (paste (root, 'combined_nuclei.csv', sep='/')) %>% data.frame()

cell_lab <- read.csv('../inst/python/clono/TF_cell.csv')
all_data %>% colony_count (cell_lab=cell_lab) -> proc_dat
write.csv (proc_dat, paste (root, 'combined_proc.csv', sep='/'))
plot_colony_count (proc_dat)

# mean cell number per colony
ggplot (proc_dat, aes (x=mean_colony, y=rel_num_colony, color=condition))+
        geom_point ()+
        TBdev::theme_TB ('dotplot', rotation=0)

# take a long time to run:
# plot_colony_count (proc_dat, y_val='mean_colony')

# ----------Total number per well----------
proc_dat %>% group_by (condition, Date) %>%
        summarise (total_num = sum(rel_num_colony)) -> by_run
by_run$celltype <- cell_lab$celltype [match (by_run$condition, toupper (cell_lab$TF))]
# filter out poor quality replicates:
by_run %>% filter (! Date %in% c('run3') ) -> by_run

cluster_level <- c('GFP', 'PITX1', 'CTNNB1', 'TFAP2C', 'HMOX1', 'ZNF750',
                   'SP6', 'MSX2', 'OVOL1', 'GRHL1', 'GCM1', 'TEAD3', 'GATA2',
                   'ZFHX3', 'GATA3', 'IRX4', 'NFE2L3', 'NR2F2', 'TFEB',
                   'TRIM28', 'MAZ', 'SSRP1', 'PCBP2', 'HMGA1')
duals <- c('GATA2_CEBPA', 'GATA2_GATA3', 'MSX2_NR2F2', 'MSX2_OVOL1', 'TEAD3_GCM1')
devtools::load_all ('..', export_all=T)
cairo_pdf (paste (save_dir, 'clonogenicity_norep_box.pdf', sep='/'), width=14, height=8)
by_run %>% filter (!Date %in% c('run6', 'run8') & !condition %in% duals ) %>%
plot_colony_count (y_val='total_num', plot_type='box', paired_test=T, 
                   new_level=cluster_level,
                   plot_pval='p.signif', display_shape=F)+
ylab ('normalized number of colonies')
dev.off()

cairo_pdf (paste (save_dir, 'clonogenicity_dual_norep_box.pdf', sep='/'), width=8, height=8)
by_run %>% filter (condition %in% c(duals, 'GATA2', 'GATA3', 'MSX2', 'OVOL1',
                                    'TEAD3', 'GCM1', 'GFP') ) %>%
plot_colony_count (y_val='total_num', plot_type='box', paired_test=T,
                   plot_pval='p.signif', display_shape=F, new_level=cluster_level)+
theme (aspect.ratio=1)
dev.off()

# save data for network plot
by_run %>% group_by (condition) %>% summarise (KD=mean (total_num)) %>%
        write.csv (paste (save_dir, 'clonogenicity_KD_efficacy.csv', sep='/'))

# ----------reproducibility----------
by_run %>% group_by (condition, Date) %>% select (!celltype) %>%
        spread ('Date', 'total_num') %>% 
        ggplot (aes (x=run3, y=run2))+
        geom_point ()+ ggrepel::geom_text_repel (aes (label=condition), size=5)+
        geom_abline (slope=1, intercept=0, color='red')+
        TBdev::theme_TB ('dotplot', rotation=0)

rmat <- reproducibility_mat (by_run)
pheatmap::pheatmap (rmat)

# ----------validation----------
true_df <- read.csv (paste (root, 'clono/proc_data/validation.csv', sep='/'))
true_df$condition <- toupper (true_df$condition)
true_df$ID <- paste (true_df$condition, true_df$trial, 'run1', sep='_')
proc_dat$ID <- paste (proc_dat$condition, proc_dat$series, proc_dat$Date, sep='_')
true_df$estimate <- proc_dat$num_colony [match (true_df$ID, proc_dat$ID)]

data (format_conf, package='TBdev')
true_df$estimate [is.na (true_df$estimate)] <- 0
rmse <- (true_df$estimate - true_df$true_num)^2 %>% mean () %>% sqrt () %>% round (2)
cairo_pdf (paste (sup_save_dir, 'clono_validate.pdf', sep='/'))
ggplot (true_df, aes (x=true_num, y=estimate))+
        geom_point ()+
        ggrepel::geom_text_repel (aes(label=condition), size=format_conf$point_fontsize)+
        geom_abline (slope=1, intercept=c(0,0)) +
        TBdev::theme_TB ('dotplot', rotation=0) +xlab ('manual') + ylab ('automatic')+
        ggtitle (paste ('RMSE=', rmse))
dev.off ()

# filtered dataframe
all_data %>% dplyr::filter (cluster != -1) %>% 
        dplyr::filter (area < 400) %>%
        dplyr::count (condition, series, Date, cluster) %>%
        dplyr::filter (n>=6 & n < 250) %>% 
        tidyr::unite (ID, c('condition', 'series', 'Date', 'cluster')) %>%
        dplyr::pull(ID) -> fil_data

all_data %>% tidyr::unite (ID, c('condition', 'series', 'Date', 'cluster'), remove=F) -> all_data
filtered <- all_data [all_data$ID %in% fil_data,]
write.csv (filtered, paste (root, 'filtered_nuclei.csv', sep='/'))
