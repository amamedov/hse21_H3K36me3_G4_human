library(ggplot2)

setwd('C:/Users/Артём/Desktop/HSE/Бионформатика/Проект/Мамедов/hse21_H3K36me3_G4_human/src/')

names <- c('Homo_Li_KPDS')

outdir <- '../images/'


for (NAME in names)
{
  bed_df <- read.delim(paste0('../data/',NAME, '.bed'), 
                       as.is = TRUE, header = FALSE)
  colnames(bed_df) <- c('chrom', 'start', 'end')
  bed_df$len <- bed_df$end - bed_df$start
  
  ggplot(bed_df) +
    aes(x = len) +
    geom_histogram() +
    ggtitle(NAME, subtitle = sprintf('Number of peaks = %s', nrow(bed_df))) +
    theme_bw()
  ggsave(paste0('len_hist.', NAME, '.pdf'), path = outdir,width = 10, dpi = 600)
}

