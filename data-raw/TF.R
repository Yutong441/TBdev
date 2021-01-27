# generate a list of transcriptional factors
# download the list from http://humantfs.ccbr.utoronto.ca/download.php
# select 1 Human TFs > Full Database > .csv
# The downloaded file is called 'DatabaseExtract_v_1.01.csv'

TF <- read.csv ('config/DatabaseExtract_v_1.01.csv')
TF <- as.character (TF$HGNC.symbol)
usethis::use_data (TF, overwrite=T)
