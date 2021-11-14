# compare in vitro cell lines
library (magrittr)
#setwd ('..')
#roxygen2::roxygenise()
#setwd ('manuscript')

devtools::load_all ('..', export_all = T)
root_dir <- '/mnt/c/Users/Yutong/Documents/bioinformatics/reproduction/'
root <- paste (root_dir, 'results/', sep='/')
merge_dir <- paste (root, 'XLYBPZ_Dylan_dir', sep='/')
save_dir <- paste (root, 'manuscript/figure5', sep='/')
all_data <- readRDS(paste (merge_dir, 'merged_blastoid.rds', sep='/') )

in_vitro <- all_data [, all_data$date=='in_vitro']
markers <- find_DE_genes (in_vitro, save_dir, group.by='broad_type', label='all_vitro', method='pairwise')
sel_cells <- c('hESC', 'hTSC-OKAE', 'bTSC', 'bTB-YANA')
#markers %>% dplyr::filter (group %in% sel_cells & compare_group %in% sel_cells) %>%
#        dplyr::select (!V1) %>% write.csv (paste (save_dir, 'DE_genes_in_vitro.csv', sep='/'))
save_dir1 <- paste (root, 'manuscript/figure1', sep='/')
vivo_mark <- find_DE_genes (all_data, save_dir1, group.by='broad_type', label='pairwise')

gtable <- all_vivo_all_vitro (markers, vivo_mark, 
                              c('hESC', 'hTSC-OKAE', 'bTSC', 'bTB-YANA'), 
                              c('ICM', 'TB', 'CTB', 'STB', 'EVT')
)

write.csv (gtable, paste (save_dir, 'vitro_vivo_DE_comp.csv', sep='/'))
markers %>% dplyr::filter (group=='hESC' & compare_group=='hTSC-OKAE') %>%
        dplyr::arrange (dplyr::desc (logFC)) %>% head ()

ptable <- p_val_all_pairs (gtable, c('hTSC-OKAE', 'bTSC', 'bTB-YANA', 'hESC'),
                           c('ICM', 'TB', 'CTB', 'STB', 'EVT'))
write.csv (ptable, paste (save_dir, 'vitro_vivo_DE_pval.csv', sep='/'))

sum_table <- signif_p_val (ptable)
write.csv (sum_table, paste (save_dir, 'vitro_vivo_DE_pval_summary.csv', sep='/'))

devtools::load_all ('..', export_all = T)
pure_vitro_DE_genes (markers, vivo_mark, 
                    c('hESC', 'hTSC-OKAE', 'bTSC', 'bTB-YANA'), 
                    c('ICM', 'TB', 'CTB', 'STB', 'EVT'), select_num=20
) -> pure_vitro_all

write.csv (pure_vitro_all, paste (save_dir, 'pure_vitro_gene_selected.csv', sep='/'))  

# ----------Over representation----------
1:nrow(pure_vitro_all) %>% as.list () %>% lapply (function (i){
                pure_vitro_all [i,'gene'] %>% as.character () %>%
                        strsplit(',') %>% unlist () %>% trimws () } ) ->gene_list
names (gene_list) <- rownames (pure_vitro_all)

library (org.Hs.eg.db)
kk <- compare_cluster_enrichment (gene_list, d, org.Hs.eg.db,
                                  enrich_area='KEGG')
display_cluster_enrichment (kk, show_graph='emap', show_num=20) + ggplot2::labs (fill='')

go <- compare_cluster_enrichment (gene_list, d, org.Hs.eg.db,
                                  enrich_area='GO')
display_cluster_enrichment (go, show_graph='emap', show_num=20) + ggplot2::labs (fill='')

# ----------GSEA----------
1:length (gene_list) %>% as.list () %>% lapply (function (i){
                group1 <- gsub ('_vs_.*$', '', names (gene_list)[i])
                group2 <- gsub ('^.*_vs_', '', names (gene_list)[i])
                markers %>% dplyr::filter (group==group1 & compare_group==group2) %>%
                        dplyr::filter (feature %in% gene_list[[i]])
          }) -> fil_markers

fil_markers <- do.call (rbind, fil_markers)

devtools::load_all ('..', export_all = T)
library (org.Hs.eg.db)
psea <- run_GSEA_pairwise (markers, org.Hs.eg.db, paste (save_dir, 'GSEA_hESC_hTSC-OKAE.csv', sep='/'),
                   group1=c('hESC', 'bTB-YANA', 'bTSC', 'hTSC-OKAE'), 
                   group2=c('hESC', 'bTB-YANA', 'bTSC', 'hTSC-OKAE'))

cairo_pdf (paste (save_dir, 'GSEA_vitro_bTSC.pdf', sep='/'), width=14, height=5)
psea %>% dplyr::filter (compare_group=='bTSC') %>%
rich_forest (markers, org.Hs.eg.db, show_num=3, show_gene_labels=2,
             shrink_ratio=0.7, nudge_x=0.2, extend_x_pos=1.2, extend_x_neg=1.0) 
dev.off ()

cairo_pdf (paste (save_dir, 'GSEA_vitro_bTB.pdf', sep='/'), width=14, height=5)
psea %>% dplyr::filter (compare_group=='bTB-YANA') %>%
rich_forest (markers, org.Hs.eg.db, show_num=3, show_gene_labels=2,
             shrink_ratio=0.7, nudge_x=0.2, extend_x_pos=1.2, extend_x_neg=1.0) 
dev.off ()
