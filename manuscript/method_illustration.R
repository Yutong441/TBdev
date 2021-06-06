# This script creates the panels for illustrating the method of analysis

library (tidyverse)
devtools::load_all ('..', export_all=T)
root_dir <- '/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/'
root <- paste (root_dir, 'results/', sep='/')
merge_dir <- paste (root, 'XLYBPZ_Dylan_dir', sep='/')
save_dir <- paste (root, 'manuscript/method', sep='/')
sup_save_dir <- paste (root, 'manuscript/figureS2', sep='/')
all_data <- get(load (paste (merge_dir, 'final_merged_tb.Robj', sep='/') ))

# ----------plain GPLVM----------
epg <- read.csv (paste (sup_save_dir, 'result/STREAM_graph.csv', sep='/'))
metadata <- all_data@meta.data 
metadata %>% filter (!broad_type %in% c('EPI', 'PE')) %>% filter (!is.na (PT1)) -> metadata
new_name <- c('TB_stem', 'EVT_branch', 'STB_branch', NA)
old_name <- c('S2,S1', 'S0,S1', 'S3,S1', 'blank')
metadata %>% slice_min (PT3, n=nrow(metadata)-3) %>%
        filter (broad_type != 'uCTB')-> plot_met

cairo_pdf (paste (save_dir, 'panelB_plain.pdf', sep='/'))
dim_red_3D (plot_met, 'PT1', 'PT2', 'PT3', 'broad_type', all_theta=50,
            all_phi=0, further_repel=T, repel_force=0.5, 
            lab_just=c(0.08, 0.02, 0.02), AP=list (edge_stroke=0.3),
            hor_just=0.1, dim_vjust=4) +theme (legend.position='none')
dev.off()

# ----------GPVLM with branch names----------
cairo_pdf (paste (save_dir, 'panelC_branch_lab.pdf', sep='/'))
plot_met %>% filter (!is.na (epil_branch)) %>%
dim_red_3D ('PT1', 'PT2', 'PT3', 'epil_branch', all_theta=50,
            all_phi=0, further_repel=F, repel_force=0.5, 
            lab_just=c(0.08, 0.02, 0.02), AP=list (edge_stroke=0.3),
            hor_just=0.1, dim_vjust=4) +theme (legend.position='none')
dev.off()

# ----------GPLVM with Trajectory line----------
data (format_conf)
epg$branch_name <- partial_relevel (new_name[match (epg$branch, old_name )], format_conf$branch_order)
cairo_pdf (paste (save_dir, 'panelD_traj.pdf', sep='/'))

plot_met %>% filter (!is.na (epil_branch)) %>%
dim_red_3D_traj ('PT1', 'PT2', 'PT3', 'epil_branch', epg, 'x', 'y',
                    'z', 'branch_name', all_theta=50, all_phi=0, further_repel=F,
                    repel_force=0.1, lab_just=c(0.08, 0.02, 0.02), magnify_text=1.3, 
                    hor_just=0.1, dim_vjust=4, label_col=NULL, AP=list(edge_stroke=0.3))+
theme (legend.position='none') + scale_color_manual (values=rep('black',3))

dev.off()

# ----------by pseudotime----------
cairo_pdf (paste (save_dir, 'panelE_PT.pdf', sep='/'))
plot_met %>% filter (!is.na (epil_branch)) %>%
dim_red_3D ('PT1', 'PT2', 'PT3', 'MGP_PT', all_theta=50,
            all_phi=0, further_repel=F, repel_force=0.5, 
            lab_just=c(0.08, 0.02, 0.02), AP=list (edge_stroke=0.3),
            hor_just=0.1, dim_vjust=4) +theme (legend.position='none')
dev.off()

# ----------swictch point----------
data_dir <- paste (root_dir, 'GPLVM/hier_GP_tf1/', sep='/')
exp_mat <- fast_read_df(paste (data_dir, 'data/STREAM_data.csv', sep='/' ))
pseudotime <- fast_read_df(paste (data_dir, 'result/infer_pt_matern.csv', sep='/' ))
off_set <- min (all_data$ori_MGP_PT, na.rm=T)
exp_mat %>% as.matrix () %>% t () %>% scale () %>% data.frame () %>% 
        add_column (pseudotime=pseudotime [,'pt_mean'] - off_set) -> exp_mat

pred_all <- data.table::fread(paste (data_dir, 'result/prediction_matern_500.csv', 
                                     sep='/' )) %>% data.frame ()
pred_all$x <- pred_all$x - min (all_data$ori_MGP_PT, na.rm=T)
pred_all$branch <- c('EVT_branch', 'STB_branch')[as.factor (pred_all$branch)]

devtools::load_all ('..', export_all=T)
gene_over_pseudotime (pred_all, exp_mat, 'GADD45G', metadata,
# This script creates the panels for illustrating the method of analysis

library (tidyverse)
devtools::load_all ('..', export_all=T)
root_dir <- '/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/'
root <- paste (root_dir, 'results/', sep='/')
merge_dir <- paste (root, 'XLYBPZ_Dylan_dir', sep='/')
save_dir <- paste (root, 'manuscript/method', sep='/')
sup_save_dir <- paste (root, 'manuscript/figureS2', sep='/')
all_data <- get(load (paste (merge_dir, 'final_merged_tb.Robj', sep='/') ))

# ----------plain GPLVM----------
epg <- read.csv (paste (sup_save_dir, 'result/STREAM_graph.csv', sep='/'))
metadata <- all_data@meta.data 
metadata %>% filter (!broad_type %in% c('EPI', 'PE')) %>% filter (!is.na (PT1)) -> metadata
new_name <- c('TB_stem', 'EVT_branch', 'STB_branch', NA)
old_name <- c('S2,S1', 'S0,S1', 'S3,S1', 'blank')
metadata %>% slice_min (PT3, n=nrow(metadata)-3) %>%
        filter (broad_type != 'uCTB')-> plot_met

cairo_pdf (paste (save_dir, 'panelB_plain.pdf', sep='/'))
dim_red_3D (plot_met, 'PT1', 'PT2', 'PT3', 'broad_type', all_theta=50,
            all_phi=0, further_repel=T, repel_force=0.5, 
            lab_just=c(0.08, 0.02, 0.02), AP=list (edge_stroke=0.3),
            hor_just=0.1, dim_vjust=4) +theme (legend.position='none')
dev.off()

# ----------GPVLM with branch names----------
cairo_pdf (paste (save_dir, 'panelC_branch_lab.pdf', sep='/'))
plot_met %>% filter (!is.na (epil_branch)) %>%
dim_red_3D ('PT1', 'PT2', 'PT3', 'epil_branch', all_theta=50,
            all_phi=0, further_repel=F, repel_force=0.5, 
            lab_just=c(0.08, 0.02, 0.02), AP=list (edge_stroke=0.3),
            hor_just=0.1, dim_vjust=4) +theme (legend.position='none')
dev.off()

# ----------GPLVM with Trajectory line----------
data (format_conf)
epg$branch_name <- partial_relevel (new_name[match (epg$branch, old_name )], format_conf$branch_order)
cairo_pdf (paste (save_dir, 'panelD_traj.pdf', sep='/'))
plot_met %>% filter (!is.na (epil_branch)) %>%
dim_red_3D_traj ('PT1', 'PT2', 'PT3', 'epil_branch', epg, 'x', 'y',
                    'z', 'branch_name', all_theta=50, all_phi=0, further_repel=F,
                    repel_force=0.1, lab_just=c(0.08, 0.02, 0.02), magnify_text=1.3, 
                    hor_just=0.1, dim_vjust=4, label_col=NULL, AP=list(edge_stroke=0.3))+
theme (legend.position='none') + scale_color_manual (values=rep('black',3))
dev.off()

# ----------by pseudotime----------
cairo_pdf (paste (save_dir, 'panelE_PT.pdf', sep='/'))
plot_met %>% filter (!is.na (epil_branch)) %>%
dim_red_3D ('PT1', 'PT2', 'PT3', 'MGP_PT', all_theta=50,
            all_phi=0, further_repel=F, repel_force=0.5, 
            lab_just=c(0.08, 0.02, 0.02), AP=list (edge_stroke=0.3),
            hor_just=0.1, dim_vjust=4) 
dev.off()

# ----------swictch point----------
data_dir <- paste (root_dir, 'GPLVM/hier_GP_tf1/', sep='/')
exp_mat <- fast_read_df(paste (data_dir, 'data/STREAM_data.csv', sep='/' ))
pseudotime <- fast_read_df(paste (data_dir, 'result/infer_pt_matern.csv', sep='/' ))
off_set <- min (all_data$ori_MGP_PT, na.rm=T)
exp_mat %>% as.matrix () %>% t () %>% scale () %>% data.frame () %>% 
        add_column (pseudotime=pseudotime [,'pt_mean'] - off_set) -> exp_mat

pred_all <- data.table::fread(paste (data_dir, 'result/prediction_matern_500.csv', 
                                     sep='/' )) %>% data.frame ()
pred_all$x <- pred_all$x - min (all_data$ori_MGP_PT, na.rm=T)
pred_all$branch <- c('EVT_branch', 'STB_branch')[as.factor (pred_all$branch)]

# peak data
save_dir2 <- paste (root, 'manuscript/figure2', sep='/')
peak_plot <- read.csv(paste (save_dir2, 'Switch_EVT.csv', sep='/'))

cairo_pdf (paste (save_dir, 'panelF_HLAG.pdf', sep='/'))
gene_over_pseudotime (pred_all, exp_mat, 'HLA.G', metadata,peak_data=peak_plot,
                      color_feature = 'broad_type', num_col=2, plot_line=T
                      )+theme (strip.text=element_blank(), legend.position='none')
dev.off()
