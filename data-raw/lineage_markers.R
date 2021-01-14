# Generate lineage markers genes
lineage_markers <- list (
        TB = c('KRT7', 'KRT18', 'TFAP2C', 'TFAP2A', 'GATA3', 'GATA2',
                'ELF5', 'CDX2', 'TEAD4', 'ETS2', 'HAND1'),
        # ? ELF5, EOMES
        EVT = c('ITGA1', 'ITGA5', 'HLA-G', 'MMP2', 'MMP9'),
        STB = c('CGA', 'TBX3', 'CSH1', 'CGB1', 'CGB2', 'LHB', 'PGF'),
        CTB = c('TP63', 'OVOL1', 'CDH1', 'CK7'),

        #EPI = c('SOX2', 'POU5F1', 'PRDM14', 'GDF3', 'TDGF1', 'KLF17', 'NODAL',
        #        'DPPA5', 'ARGFX', 'IFITM1', 'PRICKLE1', 'NANOG'),
        EPI = c('SOX2', 'POU5F1', 'NANOG', 'KLF17', 'PRDM14'),
        #PE = c('SERPINH1', 'PGFR2', 'LAMA4', 'HNF1B', 'PDGFRA', 'KIT', 'BMP2',
        #       'LBH', 'FGFR2', 'SOX17')
        PE = c( 'GATA6', 'NID2', 'HNF1B', 'PDGFRA', 'SOX17'),
        STR = c('CD14', 'CD34', 'CD90', 'CD106', 'HLA-A', 'HLA-B')
)

for (i in 1:length(lineage_markers)){
        names (lineage_markers [[i]]) <- rep (names (lineage_markers)[i],
                                              length (lineage_markers[[i]]))
}
lineage_markers <- Reduce (c, lineage_markers)
usethis::use_data (lineage_markers, overwrite=T)
