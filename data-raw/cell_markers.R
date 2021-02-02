# cell types
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

orders <- list (
        cell_order = c('Oocyte', 'Zy', '2C', '4C', '8C', 'cMor', 'Blast',
                        'cleavage', 'ICM', 'EPI', 'PSA-EPI','PE', 'HYP', 'TB',
                        'TE', 'CTB', 'STB', 'EVT', 'STR'),
        time_order = c('early', 'intermediate', 'intermediate1', 'intermediate2', 'advanced'),
        branch_order = c('TB_stem', 'STB_branch', 'EVT_branch')
)
usethis::use_data (CT, orders, overwrite=T)
