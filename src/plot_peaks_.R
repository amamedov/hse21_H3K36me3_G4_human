library(ggplot2)

#setwd('C:/Users/����/Desktop/HSE/�������������/������/�������/')

names <- c('H3K36me3_GM12878.ENCFF432EMI.hg38',
            'H3K36me3_GM12878.ENCFF475QVQ.hg38',
            'H3K36me3_GM12878.ENCFF432EMI.hg19',
            'H3K36me3_GM12878.ENCFF475QVQ.hg19')

outdir <- '../images/'


for (NAME in names)
{
  bed_df <- read.delim(paste0('../data/',NAME, '.bed'), 
                       as.is = TRUE, header = FALSE)
  colnames(bed_df) <- c('chrom', 'start', 'end', 'name', 'score')
  bed_df$len <- bed_df$end - bed_df$start
  
  ggplot(bed_df) +
    aes(x = len) +
    geom_histogram() +
    ggtitle(NAME, subtitle = sprintf('Number of peaks = %s', nrow(bed_df))) +
    theme_bw()
  ggsave(paste0('len_hist.', NAME, '.pdf'), path = outdir,width = 10, dpi = 600)
}

