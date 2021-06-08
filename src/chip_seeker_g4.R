#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
#BiocManager::install("ChIPseeker")
#BiocManager::install("clusterProfiler")
#BiocManager::install("org.Hs.eg.db")

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(clusterProfiler)
setwd('C:/Users/Артём/Desktop/HSE/Бионформатика/Проект/Ямпольская/project_x/src/')
###
data_dir <- '../data/'
out_dir <- '../images/'
name_1 <- 'Homo_Li_KPDS'
BED_FN_1 <- paste0(data_dir, name_1, '.bed')

###

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

peakAnno_1 <- annotatePeak(BED_FN_1, tssRegion=c(-3000, 3000), 
                         TxDb=txdb, annoDb="org.Hs.eg.db")

pdf(paste0(out_dir, 'chip_seeker.', name_1, '.plotAnnoPie.pdf'))
plotAnnoPie(peakAnno_1)
dev.off()
