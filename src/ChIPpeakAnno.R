#BiocManager::install("ChIPpeakAnno")
# BiocManager::install("org.Hs.eg.db")
# BiocManager::install("org.Mm.eg.db")

library(ChIPpeakAnno)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(dplyr)
#library(TxDb.Mmusculus.UCSC.mm10.knownGene)
#library(org.Mm.eg.db)

###
setwd('C:/Users/Артём/Desktop/HSE/Бионформатика/Проект/Мамедов/hse21_H3K36me3_G4_human/src/')

peaks <- toGRanges(paste0('../data/', 'H3K36me3_GM12878.intersect_with_Homo_Li_KPDS.bed')
                   , format="BED")
peaks[1:2]

annoData <- toGRanges(TxDb.Hsapiens.UCSC.hg19.knownGene)
annoData[1:2]


anno <- annotatePeakInBatch(peaks, AnnotationData=annoData, 
                            output="overlapping", 
                            FeatureLocForDistance="TSS",
                            bindingRegion=c(-2000, 300))
data.frame(anno) %>% head()

anno$symbol <- xget(anno$feature, org.Hs.egSYMBOL)
data.frame(anno) %>% head()

anno_df <- data.frame(anno)
write.table(anno_df, 
            file=paste0('../data/', 'H3K36me3_GM12878.intersect_with_Homo_Li_KPDS.genes.txt'),
            col.names = TRUE, row.names = FALSE, sep = '\t', quote = FALSE)

uniq_genes_df <- unique(anno_df['symbol'])
write.table(uniq_genes_df, 
            file=paste0('../data/', 'H3K36me3_GM12878.intersect_with_Homo_Li_KPDS.genes_uniq.txt'),
            col.names = FALSE, row.names = FALSE, sep = '\t', quote = FALSE)
