library(ggplot2)
library(dplyr)

setwd('C:/Users/Артём/Desktop/HSE/Бионформатика/Проект/Мамедов/hse21_H3K36me3_G4_human/src/')
name_1 <- 'H3K36me3_GM12878.ENCFF432EMI.hg19'
name_2 <- 'H3K36me3_GM12878.ENCFF475QVQ.hg19'
outdir <- '../images/'

## FIRST FILE FILTER
bed_df <- read.delim(paste0('../data/',name_1, '.bed'), 
                       as.is = TRUE, header = FALSE)
colnames(bed_df) <- c('chrom', 'start', 'end', 'name', 'score')
bed_df$len <- bed_df$end - bed_df$start

ggplot(bed_df) +
    aes(x = len) +
    geom_histogram() +
    ggtitle(name_1, subtitle = sprintf('Number of peaks = %s', nrow(bed_df))) +
    theme_bw()
ggsave(paste0('filter_peaks.', name_1, '.init.pdf'), path = outdir
         ,width = 10, dpi = 600)

head(bed_df %>%
  arrange(-len))

bed_df <- bed_df %>%
    arrange(-len) %>%
    filter(len < 4000)

ggplot(bed_df) +
    aes(x = len) +
    geom_histogram() +
    ggtitle(name_1, subtitle = sprintf('Number of peaks = %s', nrow(bed_df))) +
    theme_bw()
  ggsave(paste0('filter_peaks.', name_1, '.filtered.pdf'), path = outdir
         ,width = 10, dpi = 600)
  
bed_df %>%
    select(-len) %>%
    write.table(file=paste0('../data/', name_1 ,'.filtered.bed'),
                col.names = FALSE, row.names = FALSE, sep = '\t', quote = FALSE)

## SECOND FILE FILTER
bed_df <- read.delim(paste0('../data/',name_2, '.bed'), 
                     as.is = TRUE, header = FALSE)
colnames(bed_df) <- c('chrom', 'start', 'end', 'name', 'score')
bed_df$len <- bed_df$end - bed_df$start

ggplot(bed_df) +
  aes(x = len) +
  geom_histogram() +
  ggtitle(name_2, subtitle = sprintf('Number of peaks = %s', nrow(bed_df))) +
  theme_bw()
ggsave(paste0('filter_peaks.', name_2, '.init.pdf'), path = outdir
       ,width = 10, dpi = 600)

head(bed_df %>%
       arrange(-len))

bed_df <- bed_df %>%
  arrange(-len) %>%
  filter(len < 4000)

ggplot(bed_df) +
  aes(x = len) +
  geom_histogram() +
  ggtitle(name_2, subtitle = sprintf('Number of peaks = %s', nrow(bed_df))) +
  theme_bw()
ggsave(paste0('filter_peaks.', name_2, '.filtered.pdf'), path = outdir
       ,width = 10, dpi = 600)

bed_df %>%
  select(-len) %>%
  write.table(file=paste0('../data/', name_2 ,'.filtered.bed'),
              col.names = FALSE, row.names = FALSE, sep = '\t', quote = FALSE)