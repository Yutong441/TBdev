# process west 2019 data: https://pubmed-ncbi-nlm-nih-gov.ezp.lib.cam.ac.uk/31636193/

root <- "/mnt/c/users/Yutong/Documents/bioinformatics/reproduction/"
data_dir <- paste (root, 'data/West_2019', sep='')
save_robj <- 'West_R.Robj'
sc_counts <- read.table (paste (data_dir, 'featurecountsGSM3735309.txt', 
                        sep='/'), sep='\t', header=T)

# clean table
gene_name <- sc_counts$Geneid
sc_counts <- sc_counts [, -(1:6)]
rownames (sc_counts) <- gene_name

# label data
meta <- read.csv (paste (data_dir, 'Meta.csv', sep='/'), header=T)
meta$Date <- gsub ('_.*', '', meta$Meta) 
rownames (meta) <- meta$ID
dataset <- Seurat::CreateSeuratObject (sc_counts, meta.data = meta)
saveRDS (dataset, file = paste(data_dir, save_robj, sep='/'))
